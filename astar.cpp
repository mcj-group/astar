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

#include "MultiBucketQueue.h"
#include "BucketStructs.h"
#include "MultiQueue.h"
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

constexpr static uint64_t FSCORE_MASK = 0xffffffff;

using PQElement = std::tuple<uint32_t, uint32_t>;
struct stat {
  uint32_t iter = 0;
  uint32_t emptyWork = 0;
};

template<typename MQ>
void MQThreadTask(const Vertex* graph, MQ &wl, stat *stats,
                    std::atomic<uint64_t> *data, 
                    uint32_t sourceNode, uint32_t targetNode)
{
    uint32_t iter = 0UL;
    uint32_t emptyWork = 0UL;
    uint32_t gScore;
    uint32_t src;
    wl.initTID();

    while (true) {
        auto item = wl.pop();
        if (item) std::tie(gScore, src) = item.get();
        else break;

        uint64_t targetData = data[targetNode].load(std::memory_order_seq_cst);
        uint32_t targetDist = targetData & FSCORE_MASK;
        ++iter;

        // With the astar definition, our heuristic
        // will always overestimate. If the current task's
        // gScore is already greater than the targetDist,
        // then it won't ever lead to a shorter path to target.
        if (targetDist <= gScore) {
            ++emptyWork;
            continue;
        }

        uint64_t srcData = data[src].load(std::memory_order_seq_cst);
        uint32_t fScore = srcData & FSCORE_MASK;

        for (uint32_t e = 0; e < graph[src].adj.size(); e++) {
            auto& adjNode = graph[src].adj[e];
            uint32_t dst = adjNode.n;
            uint32_t nFScore = fScore + adjNode.d_cm;
            if (targetDist <= nFScore) continue;
            uint64_t dstData = data[dst].load(std::memory_order_seq_cst);

            // try CAS the neighbor with the new actual distance
            bool swapped = false;
            do {
                uint32_t dstDist = dstData & FSCORE_MASK;
                if (dstDist <= nFScore) break;
                uint64_t srcShift = src;
                uint64_t swapVal = (srcShift << 32) | nFScore;
                swapped = data[dst].compare_exchange_weak(
                    dstData, swapVal,
                    std::memory_order_seq_cst,
                    std::memory_order_seq_cst);
            } while(!swapped);
            if (!swapped) continue;

            // compute new heuristic of the neighbor
            uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[dst], &graph[targetNode]));
            if (targetDist <= nGScore) continue;

            // only push if relaxing this vertex is profitable
            wl.push(nGScore, dst);
        }
    }

    stats->iter = iter;
    stats->emptyWork = emptyWork;
}

template <typename MQ_Type>
void spawnTasks(const Vertex* graph, MQ_Type &wl, std::atomic<uint64_t> *data, 
                uint32_t threadNum, uint32_t sourceNode, uint32_t targetNode) {

    stat stats[threadNum];

    auto begin = std::chrono::high_resolution_clock::now();
    std::vector<std::thread*> workers;
    cpu_set_t cpuset;
    for (int i = 1; i < threadNum; i++) {
        CPU_ZERO(&cpuset);
        uint32_t coreID = i;
        CPU_SET(coreID, &cpuset);
        std::thread *newThread = new std::thread(
            MQThreadTask<MQ_Type>, std::ref(graph), 
            std::ref(wl), &stats[i], std::ref(data),
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
    MQThreadTask<MQ_Type>(graph, wl, &stats[0], data, sourceNode, targetNode);
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

    uint64_t totalIter = 0;
    uint64_t totalEmptyWork = 0;
    for (int i = 0; i < threadNum; i++) {
        totalIter += stats[i].iter;
        totalEmptyWork += stats[i].emptyWork;
    }
    std::cout << "empty work " << totalEmptyWork << "\n";
    std::cout << "runtime_ms " << ms << "\n";
}

void astarMQ(Vertex* graph, std::string qType, uint32_t numNodes, 
             uint32_t sourceNode, uint32_t targetNode, 
             uint32_t threadNum, uint32_t queueNum, 
             uint32_t batchSizePop, uint32_t batchSizePush,
             uint32_t delta, uint32_t bucketNum)
{
    std::atomic<uint64_t> *data = new std::atomic<uint64_t>[numNodes];
    for (uint i = 0; i < numNodes; i++) {
        data[i].store(UINT64_MAX, std::memory_order_seq_cst);
    }
    data[sourceNode] = FSCORE_MASK << 32; // source has no parent

    std::function<void(uint32_t)> prefetcher = [&] (uint32_t v) -> void {
        __builtin_prefetch(&data[v], 0, 3);
    };
    
    if (qType == "MQBucket") {
        printf("delta: %d\n", delta);
        printf("Buckets: %d\n", bucketNum);
        std::function<mbq::BucketID(uint32_t)> getBucketID = [&] (uint32_t v) -> mbq::BucketID {
            uint64_t d = data[v].load(std::memory_order_seq_cst);
            uint32_t fScore = d & FSCORE_MASK;
            uint32_t gScore = fScore + dist(&graph[v], &graph[targetNode]);
            return mbq::BucketID(gScore) >> delta;
        };
        using MQ_Bucket = mbq::MultiBucketQueue<
            decltype(getBucketID), decltype(prefetcher), 
            std::greater<mbq::BucketID>, uint32_t, uint32_t, true, true
        >;
        MQ_Bucket wl(getBucketID, prefetcher, queueNum, threadNum, delta, 
                     bucketNum, batchSizePop, batchSizePush, mbq::increasing);
        wl.push(0, sourceNode);
        spawnTasks<MQ_Bucket>(graph, wl, data, threadNum, sourceNode, targetNode);  

    } else {
        // qType == MQ
        using MQ = mbq::MultiQueue<decltype(prefetcher), std::greater<PQElement>, uint32_t, uint32_t>;
        MQ wl(prefetcher, queueNum, threadNum, batchSizePop, batchSizePush);
        wl.push(0, sourceNode);
        spawnTasks<MQ>(graph, wl, data, threadNum, sourceNode, targetNode);  
    }

    // trace back the path for verification
    uint32_t cur = targetNode;
    while (cur != sourceNode) {
        uint32_t parent = data[cur].load(std::memory_order_seq_cst) >> 32;
        graph[cur].prev = parent;
        cur = parent;
    }

    delete [] data;
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

    std::vector<uint32_t> queuesizes;

    std::vector<uint32_t> prios(numNodes, UINT32_MAX);
    prios[sourceNode] = 0;
    wl.push(std::make_tuple(0, sourceNode));

    auto begin = std::chrono::high_resolution_clock::now();

    uint32_t gScore;
    uint32_t src;
    uint32_t iter = 0;
    uint32_t emptyWork = 0;
    uint32_t maxSize = 0;
    while (!wl.empty()) {
        if (wl.size() > maxSize) maxSize = wl.size();
        queuesizes.push_back(wl.size());
        std::tie(gScore, src) = wl.top();
        wl.pop();

        uint32_t targetDist = prios[targetNode];
        uint32_t fScore = prios[src];
        ++iter;

        // With the astar definition, our heuristic
        // will always overestimate. If the current task's
        // gScore is already greater than the targetDist,
        // then it won't ever lead to a shorter path to target.
        if (targetDist <= gScore) {
            ++emptyWork;
            continue;
        }

        for (uint32_t e = 0; e < graph[src].adj.size(); e++) {
            auto& adjNode = graph[src].adj[e];
            uint32_t dst = adjNode.n;
            uint32_t nFScore = fScore + adjNode.d_cm;
            if (targetDist <= nFScore) continue;
            uint32_t d = prios[dst];
            if (d <= nFScore) continue;

            // compute new heuristic of the neighbor
            uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[dst], &graph[targetNode]));
            if (targetDist <= nGScore) continue;

            // only push if relaxing this vertex is profitable
            prios[dst] = nFScore;
            graph[dst].prev = src;
            wl.push({nGScore, dst});
        }
    }


    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
    std::cout << "runtime_ms " << ms << "\n";
    std::cout << "iter " << iter << "\n";
    std::cout << "empty work " << emptyWork << "\n";
    std::cout << "max size " << maxSize << "\n";

    uint64_t sum = 0;
    for (uint i = 0; i < queuesizes.size(); i++) {
        sum += queuesizes[i];
    }
    double s = sum;
    double n = queuesizes.size();
    double avg = s / n;
    std::cout << "avg size = " << avg << "\n";
}

int main(int argc, const char** argv) {
    if (argc < 2) {
        printf("Usage: %s <inFile> <startNode> <endNode> ", argv[0]);
        printf("[qType threadNum queueNum batchPop batchPush delta bucketNum printFull]\n");
        printf("Types: Serial / MQ / MQBucket\n");
        return -1;
    }

    Vertex* graph;
    uint32_t numNodes;
    std::tie(graph, numNodes) = LoadGraph(argv[1]);

    uint32_t sourceNode = (argc > 2) ? std::min((uint32_t)atoi(argv[2]), numNodes-1) : 1*numNodes/10;
    uint32_t targetNode = (argc > 3) ? std::min((uint32_t)atoi(argv[3]), numNodes-1) : 9*numNodes/10;
    std::string qType = (argc > 4) ? argv[4] : "Serial";
    uint32_t threadNum = (argc > 5) ? atol(argv[5]) : 1;
    uint32_t queueNum = (argc > 6) ? atol(argv[6]) : threadNum * 4;
    uint32_t batchSizePop = (argc > 7) ? atol(argv[7]) : 1;
    uint32_t batchSizePush = (argc > 8) ? atol(argv[8]) : 1;
    uint32_t delta = (argc > 9) ? atol(argv[9]) : 10;
    uint32_t bucketNum = (argc > 10) ? atol(argv[10]) : 64;
    uint32_t printFull = (argc > 11) ? atol(argv[11]) : 0;
    printf("Finding shortest path between nodes %d and %d\n", sourceNode, targetNode);
    printf("Type: %s\n", qType.c_str());
    printf("Threads: %d\n", threadNum);
    printf("Queues: %d\n", queueNum);
    printf("batchSizePop: %d\n", batchSizePop);
    printf("batchSizePush: %d\n", batchSizePush);

    if (qType == "Serial") {
        astarSerial(graph, numNodes, sourceNode, targetNode);
    } else if (qType == "MQ" || qType == "MQBucket" ) {
        astarMQ(graph, qType, numNodes, sourceNode, targetNode, threadNum, queueNum,
                batchSizePop, batchSizePush, delta, bucketNum);
    } else {
        std::cerr << "Unrecognized type: " << qType << "\n";
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
