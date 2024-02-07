#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <tuple>
#include <vector>
#include <cassert>

#include "BucketQueueTryLock.h"
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

struct Vertex;

struct Adj {
    Vertex* n;
    uint64_t d_cm;
};

struct Vertex {
    double lat, lon;  // in RADIANS
    std::vector<Adj> adj;

    // Ephemeral state (used during search)
    Vertex* prev; // nullptr if not visited (not in closed set)
};

using PQElement = std::tuple<uint64_t, uint64_t>;
using MQ_Bucket = BucketMultiQueue<std::greater<uint64_t>, uint64_t, uint64_t>;
using MQ_IO = MultiQueueIO<std::greater<PQElement>, uint64_t, uint64_t>;
struct stat {
  uint64_t iter = 0;
  uint64_t emptyWork = 0;
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
        for (uint32_t j = 0; j < n; j++) graph[i].adj[j].n = &graph[readU()];
        for (uint32_t j = 0; j < n; j++) graph[i].adj[j].d_cm = readD()*EarthRadius_cm;
        graph[i].prev = nullptr;
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

    uint64_t adjs = 0;
    for (uint32_t i = 0; i < numNodes; i++) adjs += graph[i].adj.size();
    printf("Read %d nodes, %ld adjacencies\n", numNodes, adjs);

    return std::tie(graph, numNodes);
}

// Distance function, in cm
// noinline for now to avoid penalizing the fine-grained version, which inlines
// dist and precomputes cos(dst->lat) (one can in fact precompute this across
// all tasks, since dst === target in all calls)
uint64_t dist(const Vertex* src, const Vertex* dst) __attribute__((noinline));

uint64_t dist(const Vertex* src, const Vertex* dst) {
    // Use the haversine formula to compute the great-angle radians
    double latS = std::sin(src->lat - dst->lat);
    double lonS = std::sin(src->lon - dst->lon);
    double a = latS*latS + lonS*lonS*std::cos(src->lat)*std::cos(dst->lat);
    double c = 2*std::atan2(std::sqrt(a), std::sqrt(1-a));

    uint64_t d_cm = c*EarthRadius_cm;
    return d_cm;
}


bool done __attribute__((aligned(128)));
Vertex* target __attribute__((aligned(128)));

// #if defined ASTAR_COARSE_GRAIN  // Coarse-grain version (~800 cycles/task)
// #define PLS_SINGLE_TASKFUNC visitVertex
// #define PLS_SINGLE_TASKFUNC_ARGS uint64_t, Vertex*, Vertex*

// // fScore = actual distance from the source
// // gScore = monotonic-heuristic distance to target
// static inline void visitVertex(swarm::Timestamp gScore, uint64_t fScore,
//                                Vertex* v, Vertex* parent) {
//     if (done) return;

//     if (!v->prev) {
//         v->prev = parent;
//         if (v == target) {
//             done = true;
//         } else {
//             swarm::enqueue_all<NOHINT>(v->adj.begin(), v->adj.end(),
//                              [fScore, v] (swarm::Timestamp gScore, Adj a) {
//                 // NOTE: Because heuristic is monotonic, no need to have/check closed set
//                 if (!a.n->prev) {
//                     //assert(a.d_cm >= dist(v, a.n)); // OK
//                     uint64_t nFScore = fScore + a.d_cm;
//                     // Due to limited precision, there may be off-by-ones; if
//                     // so, inherit the parent's timestamp to avoid going back
//                     // in time. This is safe, it's typically a 1-cm-off calculation
//                     uint64_t nGScore = std::max(gScore, nFScore + dist(a.n, target));
//                     swarm::enqueue(visitVertex, nGScore,
//                                  swarm::Hint::cacheLine(a.n), nFScore, a.n, v);
//                 }
//             }, gScore);
//         }
//     }
// }

// #else  // Fine-grained version, computes distances in parallel and close to the data

// static inline void queueVertex(swarm::Timestamp, uint64_t, Vertex*, Vertex*);

// static inline void visitVertex(swarm::Timestamp gScore, uint64_t fScore,
//                                Vertex* v, Vertex* parent) {
//     if (done) return;

//     if (!v->prev) {
//         v->prev = parent;
//         if (v == target) {
//             done = true;
//         } else {
//             swarm::enqueue_all<EnqFlags(NOHINT | MAYSPEC)>(
//                              v->adj.begin(), v->adj.end(),
//                              [fScore, v] (swarm::Timestamp gScore, Adj a) {
//                 uint64_t nFScore = fScore + a.d_cm;
//                 swarm::enqueue(queueVertex, gScore,
//                              {swarm::Hint::cacheLine(a.n), EnqFlags::MAYSPEC},
//                              nFScore, a.n, v);
//             }, gScore);
//         }
//     }
// }

// static inline void queueVertex(swarm::Timestamp parentGScore, uint64_t fScore,
//                                Vertex* v, Vertex* parent) {
//     if (!v->prev) {
//         uint64_t gScore = std::max(parentGScore, fScore + dist(v, target));
//         swarm::enqueue(visitVertex, gScore,
//                      // We can non-speculatively load the "done" variable. The
//                      // concerning case is when two tasks have equal timestamp,
//                      // and one of them update "done".
//                      // 1) They are visiting the same node, in which case their
//                      // hints serialize their execution, so there aren't
//                      // concurrent accesses to "done".
//                      // 2) They are visiting different nodes, but one is about
//                      // to visit the target. It is possible that the task that
//                      // doesn't visit the target reads !done. This causes more
//                      // work than necessary, but is safe. It's equivalent to
//                      // ordering that task before the final target-reaching
//                      // task.
//                      EnqFlags(SAMEHINT | MAYSPEC),
//                      fScore, v, parent);
//     }
// }

// #endif

template<bool usePrefetch=true>
void MQIOThreadTask(const Vertex* graph, MQ_IO &wl, stat *stats,
                    std::atomic<uint64_t> *prios, 
                    uint32_t sourceNode, uint32_t targetNode) {
    uint64_t iter = 0UL;
    uint64_t emptyWork = 0UL;
    uint64_t readCnt = 0UL;
    uint64_t dist;
    uint64_t src;
    PQElement* popBatch = new PQElement[batch1];
    PQElement* pushBatch = new PQElement[batch2];

    while (true) {
        uint32_t prefetchIdx;
        auto item = wl.tryPopBatch(popBatch);
        uint64_t size;
        if (item) size = item.get();
        else break;

        if (usePrefetch) {
            for (prefetchIdx = 0; prefetchIdx < size && prefetchIdx < MAX_PREFETCH; prefetchIdx++) {
                uint64_t v = std::get<1>(popBatch[prefetchIdx]);
                if (graph[v].adj.size() > MAX_PREFETCH_DEG) continue;
                __builtin_prefetch(&prios[v], 0, 3);
            }
        }


        uint32_t targetDist = prios[targetNode].load(std::memory_order_relaxed);
        uint64_t idx = 0;
        for (uint64_t i = 0; i < size; i++) {
            if (usePrefetch) {
                if (i > 0 && (i % MAX_PREFETCH == 0)) {
                    uint32_t end = prefetchIdx + MAX_PREFETCH;
                    for (; prefetchIdx < size && prefetchIdx < end; prefetchIdx++) {
                        uint64_t v = std::get<1>(popBatch[prefetchIdx]);
                        if (graph[v].adj.size() > MAX_PREFETCH_DEG) continue;
                        __builtin_prefetch(&prios[v], 0, 3);
                    }
                }
            }

            std::tie(dist, src) = popBatch[i];
#ifdef PERF
            uint32_t srcD = getPrioData(&prios[src]);
#else
            uint32_t srcD = prios[src].load(std::memory_order_relaxed);
#endif
            ++iter;
            if (srcD < dist) {
                emptyWork++;
                continue;
            }

            for (uint32_t e = 0; e < graph[src].adj.size(); e++) {
                auto& dst = graph[v].adj[e];
                uint64_t newDist = dist + dist(graph[src], dst.n);
                if (targetDist != UINT32_MAX && newDist > targetDist) continue;              
                pushBatch[idx] = {newDist, dst};
                if (usePrefetch) {
                    if (graph[dst].adj.size() < MAX_PREFETCH_DEG)
                    __builtin_prefetch(&prios[dst], 0, 3);
                }

                idx++;
                if (idx == batch2) {
                    uint64_t k = filter(prios, pushBatch, idx);
                    if (k != 0) wl.pushBatch(k, pushBatch);
                    idx = 0;
                }
            }
        }

        if (idx > 0) {
            uint64_t k = filter(prios, pushBatch, idx);
            if (k != 0) wl.pushBatch(k, pushBatch);
        }
    }

    delete [] pushBatch;
    delete [] popBatch;
    stats->iter = iter;
    stats->emptyWork = emptyWork;
    stats->readCnt = readCnt;
}

void astarMQIO(const Vertex* graph, uint32_t numNodes, 
                uint32_t sourceNode, uint32_t targetNode, 
                uint32_t threadNum, uint32_t queueNum, 
                uint32_t batchSizePop, uint32_t batchSizePush, uint32_t opt)
{
    MQ_IO wl(queueNum, threadNum, batchSizePop);
    std::atomic<uint64_t> *prios = new std::atomic<uint64_t>[numNodes];
    for (uint i = 0; i < numNodes; i++) {
        prios[i].store(UINT32_MAX, std::memory_order_relaxed);
    }
    prios[sourceNode] = 0;
    wl.push(std::make_tuple(0, sourceNode));

    stat stats[threadNum];

    auto begin = std::chrono::high_resolution_clock::now();
    std::vector<std::thread*> workers;
    cpu_set_t cpuset;
    for (int i = 1; i < threadNum; i++) {
        CPU_ZERO(&cpuset);
        uint64_t coreID = i;
        CPU_SET(coreID, &cpuset);
        if (opt == 0) {
            std::thread *newThread = new std::thread(
                MQIOThreadTask<true>, std::ref(graph), 
                std::ref(wl), &stats[i], std::ref(prios),
                sourceNode, targetNode
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
                std::ref(wl), &stats[i], std::ref(prios),
                sourceNode, targetNode
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
        MQIOThreadTask<true>(graph, wl, &stats[0], prios, sourceNode, targetNode);
    else
        MQIOThreadTask<false>(graph, wl, &stats[0], prios, sourceNode, targetNode);
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
}


uint64_t neighDist(Vertex* v, Vertex* w) {
    for (Adj a : v->adj) if (a.n == w) return a.d_cm;
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
    std::string queueType(argv[4]);
    uint32_t threadNum = atol(argv[5]);
    uint32_t queueNum = atol(argv[6]);
    uint32_t bucketNum = atol(argv[7]);
    uint32_t batchSizePop = atol(argv[8]);
    uint32_t batchSizePush = atol(argv[9]);
    uint32_t opt = atol(argv[10]);
    printf("Finding shortest path between nodes %d and %d\n", sourceNode, targetNode);
    printf("Type: %d\n", queueType);
    printf("Threads: %d\n", queueNum);
    printf("Queues: %d\n", queueNum);
    printf("Buckets: %d\n", bucketNum);
    printf("batchSizePop: %d\n", batchSizePop);
    printf("batchSizePush: %d\n", batchSizePush);
    printf("opt: %d\n", opt);

    Vertex* source = &graph[sourceNode];
    target = &graph[targetNode];
    done = false;

    if (queueType == "MQIO") {
        astarMQIO(graph, numNodes, sourceNode, targetNode, threadNum, queueNum, batchSizePop, batchSizePush, opt);
    // } else if (queueType == "MQIOPlain") {
    //     astarMQIOPlain(graph, numNodes, sourceNode, targetNode, threadNum, queueNum);
    // } else if (queueType == "MQBucket") {
    //     astarMQBucket(graph, numNodes, sourceNode, targetNode, 
    //         threadNum, queueNum, bucketNum, 
    //         batchSizePop, batchSizePush, opt);
    } else {
        std::cerr << "Unrecognized type: " << queueType << "\n";
        return 1;
    }

    // Print the resulting path
    std::vector<Vertex*> path;
    Vertex* cur = target;
    while (true) {
        path.push_back(cur);
        if (cur == source) break;
        cur = cur->prev;
        assert(cur);
    }
    std::reverse(path.begin(), path.end());

    uint64_t totalDist_cm = 0;
    for (uint32_t i = 0; i < path.size()-1; i++) {
        uint64_t curDist_cm = neighDist(path[i], path[i+1]);
        totalDist_cm += curDist_cm;
        printf("%4d: %9ld -> %9ld | %8ld.%02ld m | %8ld.%02ld m\n", i, path[i] - graph, path[i+1] - graph,
                curDist_cm / 100, curDist_cm % 100, totalDist_cm / 100, totalDist_cm % 100);
    }

    uint64_t directDist_cm = dist(source, target);
    printf("As-the-crow-flies distance:                    %8ld.%02ld m\n", directDist_cm / 100, directDist_cm % 100);

    // Save the path coordinates in binary
    FILE* outFile = fopen("path.bin", "wb");
    for (Vertex* v : path) {
        fwrite(&v->lat, sizeof(double), 1, outFile);
        fwrite(&v->lon, sizeof(double), 1, outFile);
    }
    fclose(outFile);

    return 0;
}