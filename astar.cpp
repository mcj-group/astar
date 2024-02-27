#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <tuple>
#include <vector>
#include <cassert>
#include <thread>

#include "BucketQueueTryLock.h"
#include "BucketStructs.h"
#include "MultiQueueIO.h"

/* dsm: A*-search for road maps. Conventions:
 *  - Node positions are in (lat, lon) format
 *  - All latitudes, longitudes, and distances are in radians unless specified
 *    (distances are sometimes stored in cm, with vars having suffix _cm)
 *  - Distances are great-circle distances
 *  - Adjacencies include precomputed distances to reduce the number of intermediate points in long roads
 *  - Adjacencies are symmetric (all roads are considered 2-way) and the graph must be strongly connected
 *  - Algorithm finds minimum-distance (not minimum-time) route (we'd need road speeds for minimum time,
 *    but it's pretty doable)
 */

const double EarthRadius_cm = 637100000.0;
constexpr static const unsigned MAX_PREFETCH = 64;
constexpr static const unsigned MAX_PREFETCH_DEG = 1024;
constexpr static uint64_t FSCORE_MASK = 0xffffffff;
constexpr static bucket_id UNDER_BKT = INT64_MAX - 1;

struct Vertex;

struct Adj {
    uint32_t n;
    uint32_t d_cm;
};

struct Vertex {
    double lat, lon;  // in RADIANS
    std::vector<Adj> adj;

    // Ephemeral state (used during search)
    uint32_t prev; // UINT32_MAX if not visited (not in closed set)
};

using PQElement = std::tuple<uint32_t, uint32_t>;
using BktElement = std::tuple<bucket_id, uint32_t>;
using MQ_IO = MultiQueueIO<std::greater<PQElement>, uint32_t, uint32_t>;
struct stat {
  uint32_t iter = 0;
  uint32_t emptyWork = 0;
};


std::tuple<Vertex*, uint32_t> LoadGraph(const char* file) {
    const uint32_t MAGIC_NUMBER = 0x150842A7 + 0;  // increment every time you change the file format
    std::ifstream f;
    f.open(file, std::ios::binary);
    if (!f.is_open()) {
        printf("ERROR: Could not open input file\n");
        exit(1);
    }

    auto readU = [&]() -> uint32_t {
        union U {
            uint32_t val;
            char bytes[sizeof(uint32_t)];
        };
        U u;
        f.read(u.bytes, sizeof(uint32_t));
        assert(!f.fail());
        return u.val;
    };

    auto readD = [&]() -> double {
        union U {
            double val;
            char bytes[sizeof(double)];
        };
        U u;
        f.read(u.bytes, sizeof(double));
        assert(!f.fail());
        return u.val;
    };

    uint32_t magic = readU();
    if (magic != MAGIC_NUMBER) {
        printf("ERROR: Wrong input file format (magic number %d, expected %d)\n",
                magic, MAGIC_NUMBER);
        exit(1);
    }

    uint32_t numNodes = readU();
    printf("Reading %d nodes...\n", numNodes);

    Vertex* graph = new Vertex[numNodes];
    uint32_t i = 0;
    while (i < numNodes) {
        graph[i].lat = readD();
        graph[i].lon = readD();
        uint32_t n = readU();
        graph[i].adj.resize(n);
        for (uint32_t j = 0; j < n; j++) graph[i].adj[j].n = readU();
        for (uint32_t j = 0; j < n; j++) graph[i].adj[j].d_cm = readD()*EarthRadius_cm;
        graph[i].prev = UINT32_MAX;
        i++;
    }

    f.get();
    assert(f.eof());

#if 0
    // Print graph
    for (uint32_t i = 0; i < numNodes; i++) {
        printf("%6d: %7f %7f", i, graph[i].lat, graph[i].lon);
        for (auto a: graph[i].adj) printf(" %5ld %7f", a.n-graph, a.d);
        printf("\n");
    }
#endif

    uint32_t adjs = 0;
    for (uint32_t i = 0; i < numNodes; i++) adjs += graph[i].adj.size();
    printf("Read %d nodes, %d adjacencies\n", numNodes, adjs);

    return std::tie(graph, numNodes);
}

// Distance function, in cm
// noinline for now to avoid penalizing the fine-grained version, which inlines
// dist and precomputes cos(dst->lat) (one can in fact precompute this across
// all tasks, since dst === target in all calls)
uint32_t dist(const Vertex* src, const Vertex* dst) __attribute__((noinline));

uint32_t dist(const Vertex* src, const Vertex* dst) {
    // Use the haversine formula to compute the great-angle radians
    double latS = std::sin(src->lat - dst->lat);
    double lonS = std::sin(src->lon - dst->lon);
    double a = latS*latS + lonS*lonS*std::cos(src->lat)*std::cos(dst->lat);
    double c = 2*std::atan2(std::sqrt(a), std::sqrt(1-a));

    uint32_t d_cm = c*EarthRadius_cm;
    return d_cm;
}

#ifdef PERF
uint32_t __attribute__ ((noinline)) filter(
    std::atomic<uint64_t> *datas, PQElement* pushBatch, 
    PQElement* pushBatchSrcs, uint32_t pushSize) {
#else
inline uint32_t filter(
    std::atomic<uint64_t> *datas, PQElement* pushBatch, 
    PQElement* pushBatchSrcs, uint32_t pushSize) {
#endif
    uint32_t k = 0;
    for (uint32_t i = 0; i < pushSize; i++) {
        uint32_t nGScore = std::get<0>(pushBatch[i]);
        uint32_t dst = std::get<1>(pushBatch[i]);
        uint32_t nFScore = std::get<0>(pushBatchSrcs[i]);
        uint32_t src = std::get<1>(pushBatchSrcs[i]);
#ifdef PERF
        uint64_t dstData = getPrioData(&datas[v]);
#else
        uint64_t dstData = datas[dst].load(std::memory_order_relaxed);
#endif
        bool swapped = false;
        do {
            uint32_t dstDist = dstData & FSCORE_MASK;
            if (dstDist <= nFScore) break;
            uint64_t shift = src;
            uint64_t swapVal = (shift << 32) | nFScore;
            swapped = datas[dst].compare_exchange_weak(
                dstData, swapVal,
                std::memory_order_acq_rel,
                std::memory_order_acquire);
        } while(!swapped);
        if (!swapped) continue;
        std::get<0>(pushBatch[k]) = nGScore;
        std::get<1>(pushBatch[k]) = dst;
        k++;
    }
    return k;
}

template<bool usePrefetch=true>
void MQIOThreadTask(const Vertex* graph, MQ_IO &wl, stat *stats,
                    std::atomic<uint64_t> *datas, 
                    uint32_t sourceNode, uint32_t targetNode,
                    uint32_t batchSizePop, uint32_t batchSizePush) {
    uint32_t iter = 0UL;
    uint32_t emptyWork = 0UL;
    uint32_t gScore;
    uint32_t src;
    PQElement* popBatch = new PQElement[batchSizePop];
    PQElement* pushBatch = new PQElement[batchSizePush];
    PQElement* pushBatchSrcs = new PQElement[batchSizePush];

    while (true) {
        uint32_t prefetchIdx;
        auto item = wl.tryPopBatch(popBatch);
        uint32_t size;
        if (item) size = item.get();
        else break;

        if (usePrefetch) {
            for (prefetchIdx = 0; prefetchIdx < size && prefetchIdx < MAX_PREFETCH; prefetchIdx++) {
                uint32_t v = std::get<1>(popBatch[prefetchIdx]);
                if (graph[v].adj.size() > MAX_PREFETCH_DEG) continue;
                __builtin_prefetch(&datas[v], 0, 3);
            }
        }

        uint64_t targetData = datas[targetNode].load(std::memory_order_relaxed);
        uint32_t targetDist = targetData & FSCORE_MASK;
        uint32_t idx = 0;
        for (uint32_t i = 0; i < size; i++) {
            if (usePrefetch) {
                if (i > 0 && (i % MAX_PREFETCH == 0)) {
                    uint32_t end = prefetchIdx + MAX_PREFETCH;
                    for (; prefetchIdx < size && prefetchIdx < end; prefetchIdx++) {
                        uint32_t v = std::get<1>(popBatch[prefetchIdx]);
                        if (graph[v].adj.size() > MAX_PREFETCH_DEG) continue;
                        __builtin_prefetch(&datas[v], 0, 3);
                    }
                }
            }

            std::tie(gScore, src) = popBatch[i];
#ifdef PERF
            uint64_t srcData = getPrioData(&datas[src]);
#else
            uint64_t srcData = datas[src].load(std::memory_order_relaxed);
#endif
            uint32_t fScore = srcData & FSCORE_MASK;
            ++iter;

            for (uint32_t e = 0; e < graph[src].adj.size(); e++) {
                auto& adjNode = graph[src].adj[e];
                uint32_t dst = adjNode.n;
                uint32_t nFScore = fScore + adjNode.d_cm;
                if (targetDist != UINT32_MAX && nFScore > targetDist) continue;          
                uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[dst], &graph[targetNode]));
                pushBatch[idx] = {nGScore, dst};
                pushBatchSrcs[idx] = {nFScore, src};
                if (usePrefetch) {
                    if (graph[dst].adj.size() > MAX_PREFETCH_DEG) continue;
                    __builtin_prefetch(&datas[dst], 0, 3);
                }

                idx++;
                if (idx == batchSizePush) {
                    uint32_t k = filter(datas, pushBatch, pushBatchSrcs, idx);
                    if (k != 0) wl.pushBatch(k, pushBatch);
                    idx = 0;
                }
            }
        }

        if (idx > 0) {
            uint32_t k = filter(datas, pushBatch, pushBatchSrcs, idx);
            if (k != 0) wl.pushBatch(k, pushBatch);
        }
    }

    delete [] popBatch;
    delete [] pushBatch;
    delete [] pushBatchSrcs;
    stats->iter = iter;
    stats->emptyWork = emptyWork;
}

void astarMQIO(Vertex* graph, uint32_t numNodes, 
                uint32_t sourceNode, uint32_t targetNode, 
                uint32_t threadNum, uint32_t queueNum, 
                uint32_t batchSizePop, uint32_t batchSizePush, uint32_t opt)
{
    MQ_IO wl(queueNum, threadNum, batchSizePop);
    // union of 2x32bit
    // | 63..32 | 31..0  |
    // | parent | fscore |
    std::atomic<uint64_t> *datas = new std::atomic<uint64_t>[numNodes];
    for (uint i = 0; i < numNodes; i++) {
        datas[i].store(UINT64_MAX, std::memory_order_relaxed);
    }
    datas[sourceNode] = FSCORE_MASK << 32; // source has no parent
    wl.push(std::make_tuple(0, sourceNode));

    stat stats[threadNum];

    auto begin = std::chrono::high_resolution_clock::now();
    std::vector<std::thread*> workers;
    cpu_set_t cpuset;
    for (int i = 1; i < threadNum; i++) {
        CPU_ZERO(&cpuset);
        uint32_t coreID = i;
        CPU_SET(coreID, &cpuset);
        if (opt == 0) {
            std::thread *newThread = new std::thread(
                MQIOThreadTask<true>, std::ref(graph), 
                std::ref(wl), &stats[i], std::ref(datas),
                sourceNode, targetNode,
                batchSizePop, batchSizePush
            );
            int rc = pthread_setaffinity_np(newThread->native_handle(),
                                            sizeof(cpu_set_t), &cpuset);
            if (rc != 0) {
                std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
            }
            workers.push_back(newThread);
        } else {
            std::thread *newThread = new std::thread(
                MQIOThreadTask<false>, std::ref(graph), 
                std::ref(wl), &stats[i], std::ref(datas),
                sourceNode, targetNode,
                batchSizePop, batchSizePush
            );
            int rc = pthread_setaffinity_np(newThread->native_handle(),
                                            sizeof(cpu_set_t), &cpuset);
            if (rc != 0) {
                std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
            }
            workers.push_back(newThread);
        }
    }
    CPU_ZERO(&cpuset);
    CPU_SET(0, &cpuset);
    sched_setaffinity(0, sizeof(cpuset), &cpuset);
    if (opt == 0) 
        MQIOThreadTask<true>(graph, wl, &stats[0], datas, sourceNode, targetNode, 
                                batchSizePop, batchSizePush);
    else
        MQIOThreadTask<false>(graph, wl, &stats[0], datas, sourceNode, targetNode, 
                                batchSizePop, batchSizePush);
    for (std::thread*& worker : workers) {
        worker->join();
        delete worker;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
    if (!wl.empty()) {
        std::cout << "not empty\n";
    }

    wl.stat();
    std::cout << "runtime_ms " << ms << "\n";

    // trace back the path for verification
    uint32_t cur = targetNode;
    while (cur != sourceNode) {
        uint32_t parent = datas[cur].load(std::memory_order_relaxed) >> 32;
        graph[cur].prev = parent;
        cur = parent;
    }
}

void MQIOPlainThreadTask(const Vertex* graph, MQ_IO &wl, stat *stats,
                    std::atomic<uint64_t> *datas, 
                    uint32_t sourceNode, uint32_t targetNode)
{
    uint32_t iter = 0UL;
    uint32_t emptyWork = 0UL;
    uint32_t gScore;
    uint32_t src;

    while (true) {
        auto item = wl.tryPop();
        if (item) std::tie(gScore, src) = item.get();
        else break;

        uint64_t targetData = datas[targetNode].load(std::memory_order_relaxed);
        uint32_t targetDist = targetData & FSCORE_MASK;
#ifdef PERF
        uint64_t srcData = getPrioData(&datas[src]);
#else
        uint64_t srcData = datas[src].load(std::memory_order_relaxed);
#endif
        uint32_t fScore = srcData & FSCORE_MASK;
        ++iter;

        for (uint32_t e = 0; e < graph[src].adj.size(); e++) {
            auto& adjNode = graph[src].adj[e];
            uint32_t dst = adjNode.n;
            uint32_t nFScore = fScore + adjNode.d_cm;
            if (targetDist != UINT32_MAX && nFScore > targetDist) continue;              
            uint64_t dstData = datas[dst].load(std::memory_order_relaxed);
            bool swapped = false;
            do {
                uint32_t dstDist = dstData & FSCORE_MASK;
                if (dstDist <= nFScore) break;
                uint64_t srcShift = src;
                uint64_t swapVal = (srcShift << 32) | nFScore;
                swapped = datas[dst].compare_exchange_weak(
                    dstData, swapVal,
                    std::memory_order_acq_rel,
                    std::memory_order_acquire);
            } while(!swapped);
            if (!swapped) continue;
            uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[dst], &graph[targetNode]));
            wl.push({nGScore, dst});
        }
    }

    stats->iter = iter;
    stats->emptyWork = emptyWork;
}

void astarMQIOPlain(Vertex* graph, uint32_t numNodes, 
                uint32_t sourceNode, uint32_t targetNode, 
                uint32_t threadNum, uint32_t queueNum)
{
    MQ_IO wl(queueNum, threadNum, 1);
    // union of 2x32bit
    // | 63..32 | 31..0  |
    // | parent | fscore |
    std::atomic<uint64_t> *datas = new std::atomic<uint64_t>[numNodes];
    for (uint i = 0; i < numNodes; i++) {
        datas[i].store(UINT64_MAX, std::memory_order_relaxed);
    }
    datas[sourceNode] = FSCORE_MASK << 32; // source has no parent
    wl.push(std::make_tuple(0, sourceNode));

    stat stats[threadNum];

    auto begin = std::chrono::high_resolution_clock::now();
    std::vector<std::thread*> workers;
    cpu_set_t cpuset;
    for (int i = 1; i < threadNum; i++) {
        CPU_ZERO(&cpuset);
        uint32_t coreID = i;
        CPU_SET(coreID, &cpuset);
        std::thread *newThread = new std::thread(
            MQIOPlainThreadTask, std::ref(graph), 
            std::ref(wl), &stats[i], std::ref(datas),
            sourceNode, targetNode
        );
        int rc = pthread_setaffinity_np(newThread->native_handle(),
                                        sizeof(cpu_set_t), &cpuset);
        if (rc != 0) {
            std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
        }
        workers.push_back(newThread);
    }
    CPU_ZERO(&cpuset);
    CPU_SET(0, &cpuset);
    sched_setaffinity(0, sizeof(cpuset), &cpuset);
    MQIOPlainThreadTask(graph, wl, &stats[0], datas, sourceNode, targetNode);
    for (std::thread*& worker : workers) {
        worker->join();
        delete worker;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
    if (!wl.empty()) {
        std::cout << "not empty\n";
    }

    wl.stat();
    std::cout << "runtime_ms " << ms << "\n";

    // trace back the path for verification
    uint32_t cur = targetNode;
    while (cur != sourceNode) {
        uint32_t parent = datas[cur].load(std::memory_order_relaxed) >> 32;
        graph[cur].prev = parent;
        cur = parent;
    }
}

#ifdef PERF
uint32_t __attribute__ ((noinline)) CASfilter(
    std::atomic<uint64_t> *datas, BktElement* pushBatch, 
    PQElement* pushBatchSrcs, uint32_t pushSize) {
#else
inline uint32_t CASfilter(
    std::atomic<uint64_t> *datas, BktElement* pushBatch, 
    PQElement* pushBatchSrcs, uint32_t pushSize) {
#endif
    uint32_t k = 0;
    for (uint32_t i = 0; i < pushSize; i++) {
        bucket_id bkt = std::get<0>(pushBatch[i]);
        uint32_t dst = std::get<1>(pushBatch[i]);
        uint32_t nFScore = std::get<0>(pushBatchSrcs[i]);
        uint32_t src = std::get<1>(pushBatchSrcs[i]);
#ifdef PERF
        uint64_t dstData = getPrioData(&datas[v]);
#else
        uint64_t dstData = datas[dst].load(std::memory_order_relaxed);
#endif
        bool swapped = false;
        do {
            uint32_t dstDist = dstData & FSCORE_MASK;
            if (dstDist <= nFScore) break;
            uint64_t shift = src;
            uint64_t swapVal = (shift << 32) | nFScore;
            swapped = datas[dst].compare_exchange_weak(
                dstData, swapVal,
                std::memory_order_acq_rel,
                std::memory_order_acquire);
        } while(!swapped);
        if (!swapped) continue;
        std::get<0>(pushBatch[k]) = bkt;
        std::get<1>(pushBatch[k]) = dst;
        k++;
    }
    return k;
}

template<bool usePrefetch=true, typename MQ_Bucket>
void MQBucketThreadTask(const Vertex* graph, MQ_Bucket &wl, stat *stats,
                    std::atomic<uint64_t> *datas, 
                    uint32_t sourceNode, uint32_t targetNode, uint32_t delta,
                    uint32_t batchSizePop, uint32_t batchSizePush) {
    uint32_t iter = 0UL;
    uint32_t emptyWork = 0UL;
    uint32_t gScore;
    uint32_t src;
    uint32_t* popBatch = new uint32_t[batchSizePop];
    BktElement* pushBatch = new BktElement[batchSizePush];
    PQElement* pushBatchSrcs = new PQElement[batchSizePush];

    while (true) {
        uint32_t prefetchIdx;
        auto item = wl.tryPopBatch(popBatch);
        uint32_t size;
        bucket_id poppedBkt;
        if (item) std::tie(size, poppedBkt) = item.get();
        else break;

        if (usePrefetch) {
            for (prefetchIdx = 0; prefetchIdx < size && prefetchIdx < MAX_PREFETCH; prefetchIdx++) {
                uint32_t v = popBatch[prefetchIdx];
                if (graph[v].adj.size() > MAX_PREFETCH_DEG) continue;
                __builtin_prefetch(&datas[v], 0, 3);
            }
        }

        uint64_t targetData = datas[targetNode].load(std::memory_order_relaxed);
        uint32_t targetDist = targetData & FSCORE_MASK;
        uint32_t idx = 0;
        for (uint32_t i = 0; i < size; i++) {
            if (usePrefetch) {
                if (i > 0 && (i % MAX_PREFETCH == 0)) {
                    uint32_t end = prefetchIdx + MAX_PREFETCH;
                    for (; prefetchIdx < size && prefetchIdx < end; prefetchIdx++) {
                        uint32_t v = popBatch[prefetchIdx];
                        if (graph[v].adj.size() > MAX_PREFETCH_DEG) continue;
                        __builtin_prefetch(&datas[v], 0, 3);
                    }
                }
            }

            src = popBatch[i];
#ifdef PERF
            uint64_t srcData = getPrioData(&datas[src]);
#else
            uint64_t srcData = datas[src].load(std::memory_order_relaxed);
#endif
            uint32_t fScore = srcData & FSCORE_MASK;
            uint32_t gScore = (poppedBkt == UNDER_BKT) ? fScore : (poppedBkt << delta);
            ++iter;

            for (uint32_t e = 0; e < graph[src].adj.size(); e++) {
                auto& adjNode = graph[src].adj[e];
                uint32_t dst = adjNode.n;
                uint32_t nFScore = fScore + adjNode.d_cm;
                if (targetDist != UINT32_MAX && nFScore > targetDist) continue;          
                uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[dst], &graph[targetNode]));
                pushBatch[idx] = {nGScore >> delta, dst};
                pushBatchSrcs[idx] = {nFScore, src};
                if (usePrefetch) {
                    if (graph[dst].adj.size() > MAX_PREFETCH_DEG) continue;
                    __builtin_prefetch(&datas[dst], 0, 3);
                }

                idx++;
                if (idx == batchSizePush) {
                    uint32_t k = CASfilter(datas, pushBatch, pushBatchSrcs, idx);
                    if (k != 0) wl.pushBatch(k, pushBatch);
                    idx = 0;
                }
            }
        }

        if (idx > 0) {
            uint32_t k = CASfilter(datas, pushBatch, pushBatchSrcs, idx);
            if (k != 0) wl.pushBatch(k, pushBatch);
        }
    }

    delete [] popBatch;
    delete [] pushBatch;
    delete [] pushBatchSrcs;
    stats->iter = iter;
    stats->emptyWork = emptyWork;
}

template<typename MQ_Bucket>
void MQBucketThreadTaskBase(const Vertex* graph, MQ_Bucket &wl, stat *stats,
                    std::atomic<uint64_t> *datas, 
                    uint32_t sourceNode, uint32_t targetNode, uint32_t delta) {
    uint32_t iter = 0UL;
    uint32_t emptyWork = 0UL;
    uint32_t src;
    bucket_id poppedBkt;

    while (true) {
        uint32_t prefetchIdx;
        auto item = wl.tryPopSingle();
        if (item) std::tie(src, poppedBkt) = item.get();
        else break;

        uint64_t targetData = datas[targetNode].load(std::memory_order_relaxed);
        uint32_t targetDist = targetData & FSCORE_MASK;
        uint32_t idx = 0;
#ifdef PERF
        uint64_t srcData = getPrioData(&datas[src]);
#else
        uint64_t srcData = datas[src].load(std::memory_order_relaxed);
#endif
        uint32_t fScore = srcData & FSCORE_MASK;
        uint32_t gScore = (poppedBkt == UNDER_BKT) ? fScore : (poppedBkt << delta);
        ++iter;

        for (uint32_t e = 0; e < graph[src].adj.size(); e++) {
            auto& adjNode = graph[src].adj[e];
            uint32_t dst = adjNode.n;
            uint32_t nFScore = fScore + adjNode.d_cm;
            if (targetDist != UINT32_MAX && nFScore > targetDist) continue;              
            uint64_t dstData = datas[dst].load(std::memory_order_relaxed);
            bool swapped = false;
            do {
                uint32_t dstDist = dstData & FSCORE_MASK;
                if (dstDist <= nFScore) break;
                uint64_t srcShift = src;
                uint64_t swapVal = (srcShift << 32) | nFScore;
                swapped = datas[dst].compare_exchange_weak(
                    dstData, swapVal,
                    std::memory_order_acq_rel,
                    std::memory_order_acquire);
            } while(!swapped);
            if (!swapped) continue;
            uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[dst], &graph[targetNode]));
            wl.pushSingle(nGScore >> delta, dst);
        }
    }

    stats->iter = iter;
    stats->emptyWork = emptyWork;
}

void astarMQBucket(Vertex* graph, uint32_t numNodes, 
                uint32_t sourceNode, uint32_t targetNode, 
                uint32_t threadNum, uint32_t queueNum, uint32_t bucketNum,
                uint32_t batchSizePop, uint32_t batchSizePush,
                uint32_t opt, uint32_t delta)
{
    // union of 2x32bit
    // | 63..32 | 31..0  |
    // | parent | fscore |
    std::atomic<uint64_t> *datas = new std::atomic<uint64_t>[numNodes];
    for (uint i = 0; i < numNodes; i++) {
        datas[i].store(UINT64_MAX, std::memory_order_relaxed);
    }
    datas[sourceNode] = FSCORE_MASK << 32; // source has no parent

    std::function<bucket_id(uint32_t)> getBucketID = [&] (uint32_t v) -> bucket_id {
        uint64_t d = datas[v].load(std::memory_order_acquire);
        uint32_t fScore = d & FSCORE_MASK;
        uint32_t gScore = fScore + dist(&graph[v], &graph[targetNode]);
        return bucket_id(gScore) >> delta;
    };
    using MQ_Bucket = BucketMultiQueueIO<decltype(getBucketID), std::greater<bucket_id>, uint32_t, uint32_t>;
    MQ_Bucket wl(getBucketID, queueNum, threadNum, delta, bucketNum, batchSizePop, increasing);
    wl.push(0, sourceNode);

    stat stats[threadNum];

    auto begin = std::chrono::high_resolution_clock::now();
    std::vector<std::thread*> workers;
    cpu_set_t cpuset;
    for (int i = 1; i < threadNum; i++) {
        CPU_ZERO(&cpuset);
        uint32_t coreID = i;
        CPU_SET(coreID, &cpuset);
        if (opt == 0) {
            std::thread *newThread = new std::thread(
                MQBucketThreadTask<true, MQ_Bucket>, std::ref(graph), 
                std::ref(wl), &stats[i], std::ref(datas),
                sourceNode, targetNode, delta,
                batchSizePop, batchSizePush
            );
            int rc = pthread_setaffinity_np(newThread->native_handle(),
                                            sizeof(cpu_set_t), &cpuset);
            if (rc != 0) {
                std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
            }
            workers.push_back(newThread);
        } else if (opt == 1) {
            std::thread *newThread = new std::thread(
                MQBucketThreadTask<false, MQ_Bucket>, std::ref(graph), 
                std::ref(wl), &stats[i], std::ref(datas),
                sourceNode, targetNode, delta,
                batchSizePop, batchSizePush
            );
            int rc = pthread_setaffinity_np(newThread->native_handle(),
                                            sizeof(cpu_set_t), &cpuset);
            if (rc != 0) {
                std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
            }
            workers.push_back(newThread);
        } else {
            // basic
            std::thread *newThread = new std::thread(
                MQBucketThreadTaskBase<MQ_Bucket>, std::ref(graph), 
                std::ref(wl), &stats[i], std::ref(datas),
                sourceNode, targetNode, delta
            );
            int rc = pthread_setaffinity_np(newThread->native_handle(),
                                            sizeof(cpu_set_t), &cpuset);
            if (rc != 0) {
                std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
            }
            workers.push_back(newThread);
        }
    }
    CPU_ZERO(&cpuset);
    CPU_SET(0, &cpuset);
    sched_setaffinity(0, sizeof(cpuset), &cpuset);
    if (opt == 0) {
        MQBucketThreadTask<true, MQ_Bucket>(graph, wl, &stats[0], datas, sourceNode, targetNode, delta,
                                batchSizePop, batchSizePush);
    } else if (opt == 1) {
        MQBucketThreadTask<false, MQ_Bucket>(graph, wl, &stats[0], datas, sourceNode, targetNode, delta,
                                batchSizePop, batchSizePush);
    } else {
        // basic
        MQBucketThreadTaskBase<MQ_Bucket>(graph, wl, &stats[0], datas, sourceNode, targetNode, delta);
    }
    for (std::thread*& worker : workers) {
        worker->join();
        delete worker;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
    if (!wl.empty()) {
        std::cout << "not empty\n";
    }

    wl.stat();
    std::cout << "runtime_ms " << ms << "\n";

    // trace back the path for verification
    uint32_t cur = targetNode;
    while (cur != sourceNode) {
        uint32_t parent = datas[cur].load(std::memory_order_relaxed) >> 32;
        graph[cur].prev = parent;
        cur = parent;
    }
}

void astarSerial(Vertex* graph, uint32_t numNodes, 
                uint32_t sourceNode, uint32_t targetNode)
{
    using PQ = std::priority_queue<
        PQElement, 
        std::vector<PQElement>,
        std::greater<PQElement>
    >;
    PQ wl;

    std::vector<uint32_t> prios(numNodes, UINT32_MAX);
    prios[sourceNode] = 0;
    wl.push(std::make_tuple(0, sourceNode));

    auto begin = std::chrono::high_resolution_clock::now();

    uint32_t gScore;
    uint32_t src;
    uint32_t iter = 0;
    uint32_t maxSize = 0;
    while (!wl.empty()) {
        if (wl.size() > maxSize) maxSize = wl.size();
        std::tie(gScore, src) = wl.top();
        wl.pop();

        uint32_t targetDist = prios[targetNode];
        uint32_t fScore = prios[src];
        ++iter;

        for (uint32_t e = 0; e < graph[src].adj.size(); e++) {
            auto& adjNode = graph[src].adj[e];
            uint32_t dst = adjNode.n;
            uint32_t nFScore = fScore + adjNode.d_cm;
            if (targetDist != UINT32_MAX && nFScore > targetDist) continue;              
            uint32_t d = prios[dst];
            if (d <= nFScore) continue;
            uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[dst], &graph[targetNode]));
            prios[dst] = nFScore;
            graph[dst].prev = src;
            wl.push({nGScore, dst});
        }
    }


    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
    std::cout << "runtime_ms " << ms << "\n";
    std::cout << "iter " << iter << "\n";
    std::cout << "max size " << maxSize << "\n";
}

uint32_t neighDist(Vertex* graph, uint32_t v, uint32_t w) {
    auto& vnode = graph[v];
    for (Adj a : vnode.adj) if (a.n == w) return a.d_cm;
    assert(false);  // w should be in v's adjacency list
}

int main(int argc, const char** argv) {
    if (argc < 2) {
        printf("Usage: %s <inFile> [startNode endNode] [qType] [thread queue bucketNum batchPop batchPush opt]\n", argv[0]);
        printf("Types: MQIO / MQIOPlain / MQBucket\n");
        return -1;
    }

    Vertex* graph;
    uint32_t numNodes;
    std::tie(graph, numNodes) = LoadGraph(argv[1]);

    uint32_t sourceNode = (argc > 2)? std::min((uint32_t)atoi(argv[2]), numNodes-1) : 1*numNodes/10;
    uint32_t targetNode = (argc > 3)? std::min((uint32_t)atoi(argv[3]), numNodes-1) : 9*numNodes/10;
    std::string algoType(argv[4]);
    uint32_t threadNum = atol(argv[5]);
    uint32_t queueNum = atol(argv[6]);
    uint32_t bucketNum = atol(argv[7]);
    uint32_t batchSizePop = atol(argv[8]);
    uint32_t batchSizePush = atol(argv[9]);
    uint32_t opt = atol(argv[10]);
    uint32_t delta = atol(argv[11]);
    uint32_t printFull = atol(argv[12]);
    printf("Finding shortest path between nodes %d and %d\n", sourceNode, targetNode);
    printf("Type: %s\n", algoType.c_str());
    printf("Threads: %d\n", threadNum);
    printf("Queues: %d\n", queueNum);
    printf("Buckets: %d\n", bucketNum);
    printf("batchSizePop: %d\n", batchSizePop);
    printf("batchSizePush: %d\n", batchSizePush);
    printf("opt: %d\n", opt);
    printf("delta: %d\n", delta);

    if (algoType == "MQIO") {
        astarMQIO(graph, numNodes, sourceNode, targetNode, threadNum, queueNum, batchSizePop, batchSizePush, opt);
    } else if (algoType == "MQIOPlain") {
        astarMQIOPlain(graph, numNodes, sourceNode, targetNode, threadNum, queueNum);
    } else if (algoType == "MQBucket") {
        astarMQBucket(graph, numNodes, sourceNode, targetNode, 
            threadNum, queueNum, bucketNum, 
            batchSizePop, batchSizePush, opt, delta);
    } else if (algoType == "Serial") {
        astarSerial(graph, numNodes, sourceNode, targetNode);
    } else {
        std::cerr << "Unrecognized type: " << algoType << "\n";
        return 1;
    }

    // Print the resulting path
    std::vector<uint32_t> path;
    uint32_t cur = targetNode;
    while (true) {
        path.push_back(cur);
        if (cur == sourceNode) break;
        cur = graph[cur].prev;
        // assert(cur);
    }
    std::reverse(path.begin(), path.end());

    uint32_t totalDist_cm = 0;
    for (uint32_t i = 0; i < path.size()-1; i++) {
        uint32_t curDist_cm = neighDist(graph, path[i], path[i+1]);
        totalDist_cm += curDist_cm;
        if (printFull)
            printf("%4d: %9d -> %9d | %8d.%02d m | %8d.%02d m\n", i, path[i], path[i+1],
                    curDist_cm / 100, curDist_cm % 100, totalDist_cm / 100, totalDist_cm % 100);
    }
    printf("total distance: %8d.%02d m\n", totalDist_cm / 100, totalDist_cm % 100);
            
    uint32_t directDist_cm = dist(&graph[sourceNode], &graph[targetNode]);
    printf("As-the-crow-flies distance: %8d.%02d m\n", directDist_cm / 100, directDist_cm % 100);

    // Save the path coordinates in binary
    // FILE* outFile = fopen("path.bin", "wb");
    // for (uint32_t vID : path) {
    //     Vertex* v = &graph[vID];
    //     fwrite(&v->lat, sizeof(double), 1, outFile);
    //     fwrite(&v->lon, sizeof(double), 1, outFile);
    // }
    // fclose(outFile);

    // Save the path in txt
    std::ofstream outFile("path.txt");
    outFile << "start: " << sourceNode << "\n";
    outFile << "target: " << targetNode << "\n";
    outFile << "total distance: " << totalDist_cm / 100 << "." << totalDist_cm % 100 << "m\n";
    outFile << "As-the-crow-flies distance: " << directDist_cm / 100 << "." << directDist_cm % 100 << "m\n";
    for (uint32_t vID : path) {
        outFile << vID << "\n";
    }
    outFile.close();


    return 0;
}