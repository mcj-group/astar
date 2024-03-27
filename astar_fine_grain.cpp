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
#include "util.h"

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

/* The parallel versions use 64 bit data to allow for single CAS of a 64 bit data.
 * It is a union of 2 x 32 bit that encapsulates parent and distance.
 *  | 63..32 | 31..0  |
 *  | parent | fscore |
 */

/* This is a fine-grain version of astar:
 *  - It enqueues 2 tasks: one for visiting neighbors, the other for computing the floating point
 *  - The tasks are enqueued with the same priority and same key, but the key's top bit
 *    indicates whether the task is a visitTask or a computeTask
 */

constexpr static uint64_t BOTTOM32_MASK = 0xffffffff;
constexpr static uint64_t COMPUTE_TASK_BIT = 0x80000000;
constexpr static uint64_t VERTEX_ID_MASK = ~COMPUTE_TASK_BIT;
constexpr static bucket_id UNDER_BKT = - 1;

using PQElement = std::tuple<uint32_t, uint64_t>;
using BktElement = std::tuple<bucket_id, uint64_t>;
struct stat {
  uint32_t iterVisit = 0;
  uint32_t iterCompute = 0;
  uint32_t emptyWork = 0;
};

template<typename MQ_IO>
void MQIOThreadTask(const Vertex* graph, MQ_IO &wl, stat *stats,
                        std::atomic<uint64_t> *datas, 
                        uint32_t sourceNode, uint32_t targetNode)
{
    uint32_t iterVisit = 0;
    uint32_t iterCompute = 0;
    uint32_t emptyWork = 0UL;
    uint64_t task;
    uint32_t gScore;
    wl.initTID();

    while (true) {
        auto item = wl.pop();
        if (item) std::tie(gScore, task) = item.get();
        else break;

        uint64_t targetData = datas[targetNode].load(std::memory_order_relaxed);
        uint32_t targetDist = targetData & BOTTOM32_MASK;

        // With the astar definition, our heuristic
        // will always overestimate. If the current task's
        // gScore is already greater than the targetDist,
        // then it won't ever lead to a shorter path to target.
        if (targetDist <= gScore) {
            ++emptyWork;
            continue;
        }

        uint32_t vertex = task >> 32;
        uint32_t data = task & BOTTOM32_MASK;
        
        if (data & COMPUTE_TASK_BIT) {
            // perform floating point computation and enqueue the neighbor
            ++iterCompute;
            uint32_t nFScore = data & VERTEX_ID_MASK;
            uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[vertex], &graph[targetNode]));
            if (targetDist <= nGScore) continue;
            uint64_t n = vertex;
            wl.push(nGScore, n << 32);

        } else {
            // visit each neighbor and enqueue with the actual distance
            // to each neighbor, leave the floating point heuristic 
            // computation to the compute task
            ++iterVisit;
            uint64_t srcData = datas[vertex].load(std::memory_order_relaxed);
            uint32_t fScore = srcData & BOTTOM32_MASK;

            for (uint32_t e = 0; e < graph[vertex].adj.size(); e++) {
                auto& adjNode = graph[vertex].adj[e];
                uint32_t dst = adjNode.n;
                uint32_t nFScore = fScore + adjNode.d_cm;
                if (targetDist <= nFScore) continue;              
                uint64_t dstData = datas[dst].load(std::memory_order_relaxed);
                bool swapped = false;
                do {
                    uint32_t dstDist = dstData & BOTTOM32_MASK;
                    if (dstDist <= nFScore) break;
                    uint64_t srcShift = vertex;
                    uint64_t swapVal = (srcShift << 32) | nFScore;
                    swapped = datas[dst].compare_exchange_weak(
                        dstData, swapVal,
                        std::memory_order_acq_rel,
                        std::memory_order_acquire);
                } while(!swapped);
                if (!swapped) continue;
                uint64_t n = dst;
                wl.push(gScore, (n << 32 | nFScore | COMPUTE_TASK_BIT));
            }
        }
    }

    stats->iterVisit = iterVisit;
    stats->iterCompute = iterCompute;
    stats->emptyWork = emptyWork;
}

template<bool usePlain> 
void astarMQIO(Vertex* graph, uint32_t numNodes, 
                uint32_t sourceNode, uint32_t targetNode, 
                uint32_t threadNum, uint32_t queueNum,
                uint32_t batchSizePop, uint32_t batchSizePush)
{

    std::atomic<uint64_t> *datas = new std::atomic<uint64_t>[numNodes];
    for (uint i = 0; i < numNodes; i++) {
        datas[i].store(UINT64_MAX, std::memory_order_relaxed);
    }
    datas[sourceNode] = BOTTOM32_MASK << 32; // source has no parent
    std::function<void(uint32_t)> prefetcher = [&] (uint32_t v) -> void {
        __builtin_prefetch(&datas[v], 0, 3);
    };
    using MQ_IO = MultiQueueIO<decltype(prefetcher), std::greater<PQElement>, uint32_t, uint64_t, usePlain>;
    if (usePlain) {
        batchSizePop = 1;
        batchSizePush = 1;
    }
    MQ_IO wl(prefetcher, queueNum, threadNum, batchSizePop, batchSizePush);
    uint64_t node = sourceNode;
    wl.push(dist(&graph[sourceNode], &graph[targetNode]), node << 32);

    stat stats[threadNum];

    auto begin = std::chrono::high_resolution_clock::now();
    std::vector<std::thread*> workers;
    cpu_set_t cpuset;
    for (int i = 1; i < threadNum; i++) {
        CPU_ZERO(&cpuset);
        uint32_t coreID = i;
        CPU_SET(coreID, &cpuset);
        std::thread *newThread = new std::thread(
            MQIOThreadTask<MQ_IO>, std::ref(graph), 
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
    MQIOThreadTask<MQ_IO>(graph, wl, &stats[0], datas, sourceNode, targetNode);
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


template<typename MQ_Bucket>
void MQBucketThreadTask(const Vertex* graph, MQ_Bucket &wl, stat *stats,
                        std::atomic<uint64_t> *datas, 
                        uint32_t sourceNode, uint32_t targetNode, uint32_t delta) {
    uint32_t iterVisit = 0;
    uint32_t iterCompute = 0;
    uint32_t emptyWork = 0UL;
    uint64_t task;
    bucket_id poppedBkt;
    wl.initTID();

    while (true) {
        auto item = wl.pop();
        if (item) std::tie(poppedBkt, task) = item.get();
        else break;
        
        uint64_t targetData = datas[targetNode].load(std::memory_order_relaxed);
        uint32_t targetDist = targetData & BOTTOM32_MASK;
        uint32_t vertex = task >> 32;
        uint32_t data = task & BOTTOM32_MASK;

        uint32_t gScore;
        if (poppedBkt == UNDER_BKT) {
            uint64_t d = datas[vertex].load(std::memory_order_acquire);
            uint32_t fScore = d & BOTTOM32_MASK;
            gScore = fScore + dist(&graph[vertex], &graph[targetNode]);
        } else {
            gScore = poppedBkt << delta;
        }

        // With the astar definition, our heuristic
        // will always overestimate. If the current task's
        // gScore is already greater than the targetDist,
        // then it won't ever lead to a shorter path to target.
        if (targetDist <= gScore) {
            ++emptyWork;
            continue;
        }

        if (data & COMPUTE_TASK_BIT) {
            // perform floating point computation and enqueue the neighbor
            ++iterCompute;
            uint32_t nFScore = data & VERTEX_ID_MASK;
            uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[vertex], &graph[targetNode]));
            if (targetDist <= nGScore) continue;
            uint64_t n = vertex;
            wl.push(nGScore, n << 32);

        } else {
            // visit each neighbor and enqueue with the actual distance
            // to each neighbor, leave the floating point heuristic 
            // computation to the compute task
            ++iterVisit;
            uint64_t srcData = datas[vertex].load(std::memory_order_relaxed);
            uint32_t fScore = srcData & BOTTOM32_MASK;

            for (uint32_t e = 0; e < graph[vertex].adj.size(); e++) {
                auto& adjNode = graph[vertex].adj[e];
                uint32_t dst = adjNode.n;
                uint32_t nFScore = fScore + adjNode.d_cm;
                if (targetDist <= nFScore) continue;              
                uint64_t dstData = datas[dst].load(std::memory_order_relaxed);
                bool swapped = false;
                do {
                    uint32_t dstDist = dstData & BOTTOM32_MASK;
                    if (dstDist <= nFScore) break;
                    uint64_t srcShift = vertex;
                    uint64_t swapVal = (srcShift << 32) | nFScore;
                    swapped = datas[dst].compare_exchange_weak(
                        dstData, swapVal,
                        std::memory_order_acq_rel,
                        std::memory_order_acquire);
                } while(!swapped);
                if (!swapped) continue;
                uint64_t n = dst;
                wl.push(poppedBkt, (n << 32 | nFScore | COMPUTE_TASK_BIT));

            }
        }
    }

    stats->iterVisit = iterVisit;
    stats->iterCompute = iterCompute;
    stats->emptyWork = emptyWork;
}

void astarMQBucket(Vertex* graph, uint32_t numNodes, 
                uint32_t sourceNode, uint32_t targetNode, 
                uint32_t threadNum, uint32_t queueNum, uint32_t bucketNum,
                uint32_t delta, uint32_t batchSizePop, uint32_t batchSizePush)
{
    std::atomic<uint64_t> *datas = new std::atomic<uint64_t>[numNodes];
    for (uint i = 0; i < numNodes; i++) {
        datas[i].store(UINT64_MAX, std::memory_order_relaxed);
    }
    datas[sourceNode] = BOTTOM32_MASK << 32; // source has no parent

    std::function<bucket_id(uint64_t)> getBucketID = [&] (uint64_t task) -> bucket_id {
        uint32_t v = task >> 32;
        uint64_t d = datas[v].load(std::memory_order_acquire);
        uint32_t fScore = d & BOTTOM32_MASK;
        uint32_t gScore = fScore + dist(&graph[v], &graph[targetNode]);
        return bucket_id(gScore) >> delta;
    };
    std::function<void(uint32_t)> prefetcher = [&] (uint32_t v) -> void {
        __builtin_prefetch(&datas[v], 0, 3);
    };
    using MQ_Bucket = BucketMultiQueueIO<decltype(getBucketID), decltype(prefetcher), std::greater<bucket_id>, uint32_t, uint64_t>;
    MQ_Bucket wl(getBucketID, prefetcher, queueNum, threadNum, delta, bucketNum, batchSizePop, batchSizePush, increasing);
    uint64_t node = sourceNode;
    bucket_id b = bucket_id(dist(&graph[sourceNode], &graph[targetNode])) >> delta;
    wl.push(b, node << 32);

    stat stats[threadNum];

    auto begin = std::chrono::high_resolution_clock::now();
    std::vector<std::thread*> workers;
    cpu_set_t cpuset;
    for (int i = 1; i < threadNum; i++) {
        CPU_ZERO(&cpuset);
        uint32_t coreID = i;
        CPU_SET(coreID, &cpuset);
        std::thread *newThread = new std::thread(
            MQBucketThreadTask<MQ_Bucket>, std::ref(graph), 
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
    CPU_ZERO(&cpuset);
    CPU_SET(0, &cpuset);
    sched_setaffinity(0, sizeof(cpuset), &cpuset);
    MQBucketThreadTask<MQ_Bucket>(graph, wl, &stats[0], datas, sourceNode, targetNode, delta);
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
    uint64_t node = sourceNode;
    wl.push(std::make_tuple(dist(&graph[sourceNode], &graph[targetNode]), node << 32));

    auto begin = std::chrono::high_resolution_clock::now();

    uint32_t gScore;
    uint64_t task;
    uint32_t emptyWork = 0;
    uint32_t iterVisit = 0;
    uint32_t iterCompute = 0;
    uint32_t maxSize = 0;
    while (!wl.empty()) {
        if (wl.size() > maxSize) maxSize = wl.size();
        std::tie(gScore, task) = wl.top();
        wl.pop();

        uint32_t targetDist = prios[targetNode];

        // With the astar definition, our heuristic
        // will always overestimate. If the current task's
        // gScore is already greater than the targetDist,
        // then it won't ever lead to a shorter path to target.
        if (targetDist <= gScore) {
            ++emptyWork;
            continue;
        }

        uint32_t vertex = task >> 32;
        uint32_t data = task & BOTTOM32_MASK;

        if (data & COMPUTE_TASK_BIT) {
            // perform floating point computation and enqueue the neighbor
            ++iterCompute;
            uint32_t nFScore = data & VERTEX_ID_MASK;
            uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[vertex], &graph[targetNode]));
            if (targetDist <= nGScore) continue;
            uint64_t n = vertex;
            wl.push({nGScore, n << 32});

        } else {
            // visit each neighbor and enqueue with the actual distance
            // to each neighbor, leave the floating point heuristic 
            // computation to the compute task
            ++iterVisit;
            uint32_t fScore = prios[vertex];
            for (uint32_t e = 0; e < graph[vertex].adj.size(); e++) {
                auto& adjNode = graph[vertex].adj[e];
                uint32_t dst = adjNode.n;
                uint32_t nFScore = fScore + adjNode.d_cm;
                if (targetDist <= nFScore) continue;
                uint32_t d = prios[dst];
                if (d <= nFScore) continue;
                prios[dst] = nFScore;
                graph[dst].prev = vertex;
                uint64_t n = dst;
                wl.push({gScore, (n << 32 | nFScore | COMPUTE_TASK_BIT)});
            }
        }
    }


    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
    std::cout << "runtime_ms " << ms << "\n";
    std::cout << "iter visit " << iterVisit << "\n";
    std::cout << "iter compute " << iterCompute << "\n";
    std::cout << "empty work " << emptyWork << "\n";
    std::cout << "max size " << maxSize << "\n";
}

int main(int argc, const char** argv) {
    if (argc < 2) {
        printf("Usage: %s <inFile> <startNode> <endNode> [qType threadNum bucketNum batchPop batchPush printFull]\n", argv[0]);
        printf("Types: Serial / MQIO / MQIOPlain / MQBucket\n");
        return -1;
    }

    Vertex* graph;
    uint32_t numNodes;
    std::tie(graph, numNodes) = LoadGraph(argv[1]);

    uint32_t sourceNode = (argc > 2) ? std::min((uint32_t)atoi(argv[2]), numNodes-1) : 1*numNodes/10;
    uint32_t targetNode = (argc > 3) ? std::min((uint32_t)atoi(argv[3]), numNodes-1) : 9*numNodes/10;
    std::string algoType = (argc > 4) ? argv[4] : "Serial";
    uint32_t threadNum = (argc > 5) ? atol(argv[5]) : 1;
    uint32_t bucketNum = (argc > 6) ? atol(argv[6]) : 64;
    uint32_t delta = (argc > 7) ? atol(argv[7]) : 0;
    uint32_t batchSizePop = (argc > 8) ? atol(argv[8]) : 1;
    uint32_t batchSizePush = (argc > 9) ? atol(argv[9]) : 1;
    uint32_t printFull = (argc > 10) ? atol(argv[10]) : 0;
    uint32_t queueNum = threadNum * QUEUES_PER_THREAD;
    printf("Finding shortest path between nodes %d and %d\n", sourceNode, targetNode);
    printf("Type: %s\n", algoType.c_str());
    printf("Threads: %d\n", threadNum);
    printf("Queues: %d\n", queueNum);
    printf("Buckets: %d\n", bucketNum);
    printf("delta: %d\n", delta);

    if (algoType == "MQIO") {
        astarMQIO<false>(graph, numNodes, sourceNode, targetNode, threadNum, 
            queueNum, batchSizePop, batchSizePush);
    } else if (algoType == "MQIOPlain") {
        astarMQIO<true>(graph, numNodes, sourceNode, targetNode, threadNum, 
            queueNum, batchSizePop, batchSizePush);
    } else if (algoType == "MQBucket") {
        astarMQBucket(graph, numNodes, sourceNode, targetNode, threadNum,
            queueNum, bucketNum, delta, batchSizePop, batchSizePush);
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