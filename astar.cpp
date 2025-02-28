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

constexpr static uint64_t FSCORE_MASK = 0x7fffffff; // exclude msb for vector signed compares, to prevent wraparound negatives
                                                    // for germany.bin, no vertex has (uint32_t dstData)[31] (MSB) set, so this is okay

/* Vector stuff ================================== */
#include <immintrin.h>

/* Knobs */
/*     Mask Option (default orig) */
// #define AWU_MASK_SCALAR
#define AWU_MASK_VECTOR

/*     Vectorization Scheme (no default) */
#define AWU_VECSCHEME_ALWAYS   // vecld: always use vector instructions to construct mask. Use inner mask to not-gather the excess vector elements
// #define AWU_VECSCHEME_FILLVEC  // dynvec angus: vectorize only if (vertex.remaining_edges >= Bee), fill Bee-wide vector with Bee/(8) loops of vgather, else scalar
// #define AWU_VECSCHEME_MOSTOFVEC   // dynvec mcj: Bee-bit-wide mask, vectorize (vertex.remaining_edges // (8) ) times, scalar (vertex.remaining_edges % (8) times)

/*     Mask Width (default 32) (does not affect AWU_VECSCHEME_ALWAYS) */
// #define AWU_MASKWIDTH_64

#ifdef AWU_MASKWIDTH_64
constexpr uint32_t B = 64;   // 8 avx2 instructions long: avx2: 256b / uint32_t = (8)
constexpr __uint128_t UINT128_MAX =__uint128_t(__int128_t(-1L));
#elif defined(AWU_VECSCHEME_ALWAYS)
constexpr uint32_t B = 8;   // 1 avx2 instructions long
#else
constexpr uint32_t B = 32;   // 4 avx2 instructions long
#endif

/* Functions */
static inline __m256i shiftToMask256(uint8_t shift) {
  // for the lower <shift> elements, set most sig bit to 1
  int mask = 0xffu >> shift;
  // shift each bit of the mask into the msb of its corresponding element
  __m256i vmask = _mm256_set1_epi32(mask);                 // copy mask into all 32b elements
  const __m256i shift_count = _mm256_setr_epi32(31, 30, 29, 28, 27, 26, 25, 24); // reverse order: MSB [24,...,31] LSB
  __m256i msb_mask = _mm256_sllv_epi32(vmask, shift_count); // shift each mask bit up to its element's msb
  return msb_mask;
}

#ifdef AWU_MASKWIDTH_64
constexpr uint32_t countr_zero(uint64_t x) noexcept {
  return (x) ? ((uint32_t)(__builtin_ctzl(x))) : ((uint32_t)(sizeof(uint64_t) * 8));
}
#else
constexpr uint32_t countr_zero(uint32_t x) noexcept {
  return (x) ? ((uint32_t)(__builtin_ctz(x))) : ((uint32_t)(sizeof(uint32_t) * 8));
}
#endif

#include <sstream>
#include <iomanip>
std::string m256i_to_string(const char* name, __m256i vec, bool printhex=false) {
    alignas(32) int values[8];
    _mm256_store_si256(reinterpret_cast<__m256i*>(values), vec);

    std::ostringstream oss;
    oss << name << ": ";
    for (int i = 0; i < 8; ++i) {
        if (printhex) {
            oss << i << "[0x" << std::setfill('0') << std::setw(8) << std::hex << values[i] << std::dec << "] ";
        } else {
            oss << i << "[" << values[i] << "] ";
        }
    }

    // Return the resulting string
    return oss.str();
}

/* ================================================ */


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

        uint64_t targetData = data[targetNode].load(std::memory_order_relaxed);
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

        uint64_t srcData = data[src].load(std::memory_order_relaxed);
        uint32_t fScore = srcData & FSCORE_MASK;
        // std::cout << "src " << src << " =========" << std::endl;

#ifdef AWU_MASK_SCALAR // mask opt
        uint32_t eBegin = 0;
        uint32_t eEnd = graph[src].adj.size();
        constexpr uint32_t B = 8;   // 4 avx2 instructions long
        Adj * adjbase = const_cast<Adj *>(reinterpret_cast<const Adj*>(&graph[src].adj[0]));
        static_assert(sizeof(graph[src].adj[0]) == sizeof(uint64_t));

        // std::cout << "src[" << src << "] ==========" << std::endl;

        for (uint32_t e = 0; e < eEnd; e += B) {
            uint32_t mask = 0;
            for (auto f = e; f != std::min(e + B, eEnd); f++) {
                auto& adjNode = *(adjbase + f);
                uint32_t nFScore = fScore + adjNode.d_cm;
                bool cmp1 = targetDist > nFScore; // if often dont continue, then high overhead

                uint32_t dst = adjNode.n;
                uint64_t dstData = data[dst].load(std::memory_order_relaxed); // make sure vec gather preddicate on cmp1
                uint32_t dstDist = dstData & FSCORE_MASK;
                bool cmp2 = dstDist > nFScore;
                mask |= (cmp1 & cmp2) << (f - e);
            }
            // std::cout << "    mask[0x" << std::setfill('0') << std::setw(8) << std::hex << mask << std::dec << "]" << std::endl;

            for (uint32_t i = countr_zero(mask);
               i < B;
               i = countr_zero(mask & (UINT64_MAX << (i + 1)))) {
                auto& adjNode = *(adjbase + e + i);
                uint32_t nFScore = fScore + adjNode.d_cm;
                uint32_t dst = adjNode.n;
                uint64_t dstData = data[dst].load(std::memory_order_relaxed);

                // try CAS the neighbor with the new actual distance
                bool swapped = false;
                do {
                    uint32_t dstDist = dstData & FSCORE_MASK;
                    if (dstDist <= nFScore) break;
                    uint64_t srcShift = src;
                    uint64_t swapVal = (srcShift << 32) | nFScore;
                    swapped = data[dst].compare_exchange_weak(
                        dstData, swapVal,
                        std::memory_order_acq_rel,
                        std::memory_order_acquire);
                } while(!swapped);
                if (!swapped) continue;

                // compute new heuristic of the neighbor
                uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[dst], &graph[targetNode])); // delay vectorizing dist
                if (targetDist <= nGScore) continue;

                // only push if relaxing this vertex is profitable
                wl.push(nGScore, dst);
                // std::cout << "        i[" << i << "] pushed dst[" << dst << "]" << std::endl;
            }
        }
        // std::cout << std::endl;

#elif defined(AWU_MASK_VECTOR) // mask opt
        uint32_t eBegin = 0;
        uint32_t eEnd = graph[src].adj.size();
        constexpr uint32_t B = 8;   // avx2 256b vectors fit 8x 32b values

        const int *adjbase_avx = reinterpret_cast<const int*>(&graph[src].adj[0]);
        Adj *adjbase_scal = const_cast<Adj *>(reinterpret_cast<const Adj*>(&graph[src].adj[0]));
        std::atomic<uint32_t> *database = reinterpret_cast<std::atomic<uint32_t> *>(data);
        const int *database_avx = reinterpret_cast<const int*>(database);
        // *(note)
        // treat type Adj (64b struct) as {uint32_t, uint32_t}
        // auto& adjNode = *(adjbase + f);
        // type Adj is struct {uint32_t n , uint32_t d_cm} @ word 0 and word 1
        // so Adj[i] := uint64_t[i] = uint32_t[2*i] 
        static_assert(sizeof(graph[src].adj[0]) == sizeof(uint64_t));

        __m256i targetDists = _mm256_set1_epi32(targetDist);
        static_assert(sizeof(targetDist) == sizeof(uint32_t));
        __m256i fScores = _mm256_set1_epi32(fScore);
        static_assert(sizeof(fScore) == sizeof(uint32_t));

        // type Adj is struct {uint32_t n , uint32_t d_cm} @ word 0 and word 1
        const __m256i FSCORE_VMASK = _mm256_set1_epi32(FSCORE_MASK);
        const __m256i dcm_offsets = _mm256_set_epi32(15, 13, 11, 9, 7, 5, 3, 1);    // odd word indices
        const __m256i n_offsets = _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0);      // even word indices

        // std::cout << "src[" << src << "] ==========" << std::endl;

        for (uint32_t e = 0; e < eEnd; e += B) {
    #ifdef AWU_VECSCHEME_ALWAYS // vec scheme
            // type Adj is struct {uint32_t n , uint32_t d_cm} @ word 0 and word 1
            // gather a vector of d_cm's, which is every odd indexed uint32_t

            // the gist
            // mask_i64gather adj_L = Adj[i+0 : i+3] as vec<i64>
            // mask_i64gather adj_H = Adj[i+4 : i+7] as vec<i64>
            // cast adj_L and adj_H as vec<i32>
            // collect ns = {L0, L2, L4, L6, H0, H2, H4, H6}
            // collect dcms = {L1, L3, L5, L7, H1, H3, H5, H7}

            __m256i es = _mm256_set1_epi32(e*2);  // *(note), each element is 2*uint32_t
            __m256i dcm_indices = _mm256_add_epi32(es, dcm_offsets);
            // std::cout << m256i_to_string("dcm_indices", dcm_indices) << std::endl;

            uint8_t shift = ((e + B < eEnd) ? 0 : ((e+B) - eEnd));
            __m256i innermask = shiftToMask256(shift);  // uses the msb of each i32 element as maskbit

            __m256i dcms = _mm256_mask_i32gather_epi32(_mm256_setzero_si256(), adjbase_avx, dcm_offsets, innermask, 4);
            __m256i nFScores = _mm256_add_epi32(fScores, dcms);
            __m256i cmp1s = _mm256_cmpgt_epi32(targetDists, nFScores); // each arg: 1s if new < old, *compares signed ints
            cmp1s = _mm256_and_si256(cmp1s, innermask); // only care about msb<i32>, set 0 if element is unused
            // std::cout << "    " << m256i_to_string("dcms", dcms) << std::endl;
            // std::cout << "    " << m256i_to_string("targetDists", targetDists) << std::endl;
            // std::cout << "    " << m256i_to_string("nFScores", nFScores) << std::endl;
            // std::cout << "    " << m256i_to_string("cmp1s(msb)", cmp1s, true) << std::endl;
            
            __m256i n_indices = _mm256_add_epi32(es, n_offsets);
            __m256i dsts = _mm256_mask_i32gather_epi32(_mm256_setzero_si256(), adjbase_avx, n_indices, cmp1s, 4);
            // std::cout << "    " << m256i_to_string("n_indices", n_indices) << std::endl;
            // std::cout << "    " << m256i_to_string("dsts", dsts) << std::endl;

            // uint64_t dstData from accessing std::atomic<uint64_t> *data
            // then uint32_t dstDist = dstData & FSCORE_MASK; (only the lower 31 bits (int excl signed bit))
            // (1) recast data at atomic<32 bit> to get more concurrent accesses
            // (2) mask again to remove sign bit
            // (3) vector cmp (dstDist > nFScore) then make mask
            __m256i dstDatas = _mm256_mask_i32gather_epi32(_mm256_setzero_si256(), database_avx, dsts, innermask, 4);
            __m256i dstDists = _mm256_and_si256(dstDatas, FSCORE_VMASK);
            __m256i cmp2s = _mm256_cmpgt_epi32(dstDists, nFScores); // *compares signed ints
            // std::cout << "    " << m256i_to_string("dstDists", dstDists) << std::endl;
            // std::cout << "    " << m256i_to_string("cmp2s", cmp2s, true) << std::endl;
            // HERE TODO

            __m256i cmp = _mm256_and_si256(cmp1s, cmp2s);
            // std::cout << "    " << m256i_to_string("cmp", cmp, true) << std::endl;
            // std::cout << "\n" << std::endl;
            __m256 cmp_cast = _mm256_castsi256_ps(cmp);           // concat msb of each 32b int
            uint32_t mask_prep = static_cast<uint32_t>(_mm256_movemask_ps(cmp_cast));
            uint32_t mask = mask_prep & static_cast<uint32_t>(0xffu >> shift);

    #elif defined(AWU_VECSCHEME_FILLVEC) // vec scheme
        #ifdef AWU_MASKWIDTH_64
            uint64_t mask = 0;
        #else
            uint32_t mask = 0;
        #endif
            if (e+B < eEnd) {
                for (auto f = e; f < e+B; f += 8) {
                    // type Adj is struct {uint32_t n , uint32_t d_cm} @ word 0 and word 1
                    // gather a vector of d_cm's, which is every odd indexed uint32_t
                    __m256i fs = _mm256_set1_epi32(f);  // *(note)
                    __m256i offsets = _mm256_set_epi32(15, 13, 11, 9, 7, 5, 3, 1);
                    __m256i indices = _mm256_add_epi32(fs, offsets);                // offset each f by 1,3,5... to get f.d_cm

                    __m256i dcms = _mm256_mask_i32gather_epi32(_mm256_setzero_si256(), adjbase, indices, _mm256_set1_epi64x(-1), 4);

                    __m256i nFScores = _mm256_add_epi32(fScores, dcms);

                    __m256i cmp = _mm256_cmpgt_epi32(targetDists, nFScores); // each arg: 1s if new < old, *compares signed ints
                    __m256 cmp_cast = _mm256_castsi256_ps(cmp);           // concat msb of each 32b int
                    uint8_t shift = ((e + B < eEnd) ? 0 : ((e+B) - eEnd));
        #ifdef AWU_MASKWIDTH_64
                    uint64_t mask_prep = static_cast<uint64_t>(_mm256_movemask_ps(cmp_cast));
                    mask |= static_cast<uint64_t>((mask_prep & static_cast<uint64_t>(0xffu >> shift)) << (f - e));
        #else
                    uint32_t mask_prep = static_cast<uint32_t>(_mm256_movemask_ps(cmp_cast));
                    mask |= static_cast<uint32_t>((mask_prep & static_cast<uint32_t>(0xffu >> shift)) << (f - e));
        #endif
                }
            } else {
                for (auto f = e; f != std::min(e + B, eEnd); f++) {
                    const Adj & adjNode = *reinterpret_cast<const Adj*>(adjbase + f);
                    uint32_t nFScore = fScore + adjNode.d_cm;
        #ifdef AWU_MASKWIDTH_64
                    uint64_t cmp = static_cast<uint64_t>(targetDist > nFScore);
                    mask |= static_cast<uint64_t>(cmp << (f - e));
        #else
                    uint32_t cmp = static_cast<uint32_t>(targetDist > nFScore);
                    mask |= static_cast<uint32_t>(cmp << (f - e));
        #endif
                }
            }

    #elif defined(AWU_VECSCHEME_MOSTOFVEC) // vec scheme
        #ifdef AWU_MASKWIDTH_64
            uint64_t mask = 0;
        #else
            uint32_t mask = 0;
        #endif
            auto f = e;
            for(; f + 8 < std::min(e + B, eEnd); f += 8) {
                // type Adj is struct {uint32_t n , uint32_t d_cm} @ word 0 and word 1
                // gather a vector of d_cm's, which is every odd indexed uint32_t
                __m256i fs = _mm256_set1_epi32(f);
                __m256i offsets = _mm256_set_epi32(15, 13, 11, 9, 7, 5, 3, 1);
                __m256i indices = _mm256_add_epi32(fs, offsets);                // offset each f by 1,3,5... to get f.d_cm

                __m256i dcms = _mm256_mask_i32gather_epi32(_mm256_setzero_si256(), adjbase, indices, _mm256_set1_epi64x(-1), 4);

                __m256i nFScores = _mm256_add_epi32(fScores, dcms);

                __m256i cmp = _mm256_cmpgt_epi32(targetDists, nFScores); // each arg: 1s if new < old, *compares signed ints
                __m256 cmp_cast = _mm256_castsi256_ps(cmp);           // concat msb of each 32b int
                uint8_t shift = ((e + B < eEnd) ? 0 : ((e+B) - eEnd));
        #ifdef AWU_MASKWIDTH_64
                uint64_t mask_prep = static_cast<uint64_t>(_mm256_movemask_ps(cmp_cast));
                mask |= static_cast<uint64_t>((mask_prep & static_cast<uint64_t>(0xffu >> shift)) << (f - e));
        #else
                uint32_t mask_prep = static_cast<uint32_t>(_mm256_movemask_ps(cmp_cast));
                mask |= static_cast<uint32_t>((mask_prep & static_cast<uint32_t>(0xffu >> shift)) << (f - e));
        #endif
            }
            for (; f != std::min(e + B, eEnd); f++) {
                const Adj & adjNode = *reinterpret_cast<const Adj*>(adjbase + f);
                uint32_t nFScore = fScore + adjNode.d_cm;
        #ifdef AWU_MASKWIDTH_64
                uint64_t cmp = static_cast<uint64_t>(targetDist > nFScore);
                mask |= static_cast<uint64_t>(cmp << (f - e));
        #else
                uint32_t cmp = static_cast<uint32_t>(targetDist > nFScore);
                mask |= static_cast<uint32_t>(cmp << (f - e));
        #endif
            }

    #endif // vec scheme
            // Adj * adjbase = const_cast<Adj *>(reinterpret_cast<const Adj*>(&graph[src].adj[0]));
            // std::cout << "    mask[0x" << std::setfill('0') << std::setw(8) << std::hex << mask << std::dec << "]\n" << std::endl;

            // this loop is identical to the scalar one, so it is functionally correct
            for (uint32_t i = countr_zero(mask);
               i < B;
               i = countr_zero(mask & (UINT64_MAX << (i + 1)))) {
                // const Adj & adjNode = *reinterpret_cast<const Adj*>(adjbase_scal + e + i);
                auto& adjNode = *(adjbase_scal + e + i);
                uint32_t dst = adjNode.n;

                // if (targetDist <= nFScore) continue; // if does not change, can remove this cond
                uint64_t dstData = data[dst].load(std::memory_order_relaxed); // moved up from below computing nFScore

                uint32_t nFScore = fScore + adjNode.d_cm; // does this change from prev loop?
                                                          // no change in fScore or d_cm

                // try CAS the neighbor with the new actual distance
                bool swapped = false;
                do {
                    uint32_t dstDist = dstData & FSCORE_MASK;
                    if (dstDist <= nFScore) break;          // this should never be taken, since it was checked above
                    uint64_t srcShift = src;
                    uint64_t swapVal = (srcShift << 32) | nFScore;
                    swapped = data[dst].compare_exchange_weak(
                        dstData, swapVal,
                        std::memory_order_acq_rel,
                        std::memory_order_acquire);
                } while(!swapped);
                if (!swapped) continue;

                // compute new heuristic of the neighbor
                uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[dst], &graph[targetNode]));
                if (targetDist <= nGScore) continue;

                // only push if relaxing this vertex is profitable
                wl.push(nGScore, dst);
                // std::cout << "        i[" << i << "] pushed dst[" << dst << "]" << std::endl;
            }
        }
        // std::cout << std::endl;

#else
        for (uint32_t e = 0; e < graph[src].adj.size(); e++) {
            auto& adjNode = graph[src].adj[e];
            uint32_t dst = adjNode.n;
            uint32_t nFScore = fScore + adjNode.d_cm;
            if (targetDist <= nFScore) continue;  // dont check this in vector code (this is not the mask, can & this with the shift mask)
            uint64_t dstData = data[dst].load(std::memory_order_relaxed); // this is the irreg load, which I am gathering

            // try CAS the neighbor with the new actual distance
            bool swapped = false;
            do {
                uint32_t dstDist = dstData & FSCORE_MASK;
                if (dstDist <= nFScore) break; // check this in vec, double check that this is the mask I am vectorizing
                uint64_t srcShift = src;
                uint64_t swapVal = (srcShift << 32) | nFScore;
                swapped = data[dst].compare_exchange_weak(
                    dstData, swapVal,
                    std::memory_order_acq_rel,
                    std::memory_order_acquire);
            } while(!swapped);
            if (!swapped) continue;

            // compute new heuristic of the neighbor
            uint32_t nGScore = std::max(gScore, nFScore + dist(&graph[dst], &graph[targetNode])); // this line calls sinusoid (mcj thinks), expensive
            if (targetDist <= nGScore) continue;

            // only push if relaxing this vertex is profitable
            wl.push(nGScore, dst);
        }
#endif // mask opt
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
        data[i].store(UINT64_MAX, std::memory_order_relaxed);
    }
    data[sourceNode] = FSCORE_MASK << 32; // source has no parent

    std::function<void(uint32_t)> prefetcher = [&] (uint32_t v) -> void {
        __builtin_prefetch(&data[v], 0, 3);
    };
    
    if (qType == "MQBucket") {
        printf("delta: %d\n", delta);
        printf("Buckets: %d\n", bucketNum);
        std::function<mbq::BucketID(uint32_t)> getBucketID = [&] (uint32_t v) -> mbq::BucketID {
            uint64_t d = data[v].load(std::memory_order_acquire);
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
    // std::cout << "pre-verif targetNode = " << cur << std::endl;
    while (cur != sourceNode) {
        uint32_t parent = data[cur].load(std::memory_order_relaxed) >> 32;
        // std::cout << "verif parent = " << parent << std::endl;
        graph[cur].prev = parent;
        cur = parent;
    }
    // std::cout << "after verif\n" << std::endl;

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

    std::cout << "Vectorization:";
#ifdef AWU_MASK_SCALAR
    std::cout << " Scalar";
#elif defined(AWU_MASK_VECTOR)
    std::cout << " Vector";
#   ifdef AWU_VECSCHEME_ALWAYS
    std::cout << " Always";
#   elif defined(AWU_VECSCHEME_FILLVEC)
    std::cout << " FillVec";
#   elif defined(AWU_VECSCHEME_MOSTOFVEC)
    std::cout << " MostOfVec";
#   else
    std::cout << "\nError no vec scheme";
#   endif
#   ifdef AWU_MASKWIDTH_64
    std::cout << " MaskWidth-64b";
#   else
    std::cout << " MaskWidth-Default32b";
#   endif
#else
    std::cout << " Orig";
#endif
    std::cout << std::endl;

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
