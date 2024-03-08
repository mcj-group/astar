#ifndef ASTAR_UTIL
#define ASTAR_UTIL

const double EarthRadius_cm = 637100000.0;
const uint32_t QUEUES_PER_THREAD = 4;

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

uint32_t neighDist(Vertex* graph, uint32_t v, uint32_t w) {
    auto& vnode = graph[v];
    for (Adj a : vnode.adj) if (a.n == w) return a.d_cm;
    assert(false);  // w should be in v's adjacency list
}

#endif // ASTAR_UTIL
