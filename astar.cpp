#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <tuple>
#include <vector>
#include <cassert>

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

uint64_t neighDist(Vertex* v, Vertex* w) {
    for (Adj a : v->adj) if (a.n == w) return a.d_cm;
    assert(false);  // w should be in v's adjacency list
}

int main(int argc, const char** argv) {
    if (argc < 2 || argc > 4) {
        printf("Usage: %s <inFile> [startNode endNode]\n", argv[0]);
        return -1;
    }

    Vertex* graph;
    uint32_t numNodes;
    std::tie(graph, numNodes) = LoadGraph(argv[1]);

    uint32_t sourceNode = (argc > 2)? std::min((uint32_t)atoi(argv[2]), numNodes-1) : 1*numNodes/10;
    uint32_t targetNode = (argc > 3)? std::min((uint32_t)atoi(argv[3]), numNodes-1) : 9*numNodes/10;
    printf("Finding shortest path between nodes %d and %d\n", sourceNode, targetNode);

    Vertex* source = &graph[sourceNode];
    target = &graph[targetNode];
    done = false;

    // swarm::enqueue(visitVertex, dist(source, target),
    //         {swarm::Hint::cacheLine(source), EnqFlags::MAYSPEC}, 0ul, source,
    //         (Vertex*)-1ul /*can't be null, lest we revisit the source*/);
    // swarm::run();

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