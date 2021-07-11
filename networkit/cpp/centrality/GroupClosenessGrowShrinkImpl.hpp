/*
 *  GroupClosenessGrowShrinkImpl.hpp
 *
 *  Created on: 19.12.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_CENTRALITY_GROUP_CLOSENESS_GROW_SHRINK_IMPL_HPP_
#define NETWORKIT_CENTRALITY_GROUP_CLOSENESS_GROW_SHRINK_IMPL_HPP_

#ifdef __AVX2__
#include <immintrin.h>
#include <networkit/auxiliary/AlignedAllocator.hpp>
#endif // __AVX2__

#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <random>
#include <unordered_map>
#include <vector>

#include <networkit/distance/Diameter.hpp>
#include <networkit/graph/Graph.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

namespace NetworKit {
namespace GroupClosenessGrowShrinkDetails {

template <typename WeightType>
class GroupClosenessGrowShrinkImpl final {

    static constexpr WeightType infdistance = std::numeric_limits<WeightType>::max();
    // Number of packed repetitions done at once to estimate the DAG size
    static constexpr count K = 16;
    static constexpr float maxInt16 = static_cast<float>(std::numeric_limits<uint16_t>::max());

public:
    GroupClosenessGrowShrinkImpl(const Graph &G, std::vector<node> group, bool extended = false,
                                 count insertions = 0, count maxIterations = 100);

    void run();

    std::vector<node> groupMaxCloseness() const;

    count numberOfIterations() const;

private:
    const Graph *G;
    std::vector<node> group;
    // Whether or not to use the extended version of the grow-shrink algorithm.
    const bool extended;
    // Number of consecutive node insertions and removals per iteration.
    count insertions;
    // Maximum number of iterations allowed.
    const count maxIterations;

    std::vector<std::reference_wrapper<std::mt19937_64>> urngs;
    std::vector<WeightType> distance, distance_;
    std::vector<bool> visited;
    std::vector<uint32_t> sumOfMins;
    std::unordered_map<node, index> idxMap;

    std::vector<WeightType> increment;
    std::vector<node> nearest, nearest_, stack;
    std::queue<index> nextIdx;

    count totalSwaps, stackSize;

    static count computeConsecutiveInsertions(const Graph &G, size_t groupSize);

    void init();

#ifdef __AVX2__
    union RandItem {
        uint16_t items[K];
        __m256i vec;
    };
    std::vector<RandItem, AlignedAllocator<RandItem, sizeof(RandItem)>> randVec;
#else
    std::vector<uint16_t> randVec;
#endif // __AVX2__

    std::vector<std::uniform_int_distribution<uint32_t>> intDistributions;
    void initRandomVec();

    struct CompareDistance {
        CompareDistance(const std::vector<WeightType> &distance) : distance(&distance) {}

        bool operator()(const node x, const node y) const noexcept {
            return (*distance)[x] < (*distance)[y];
        }

    private:
        const std::vector<WeightType> *distance;
    };

    tlx::d_ary_addressable_int_heap<node, 2, CompareDistance> heap;
    tlx::d_ary_addressable_int_heap<node, 2, CompareDistance> heap_;

    void dijkstra();
    bool findAndSwap();
    void bfsFromGroup();
    void computeFarnessIncrement();
    // Computes real decrement of farness
    WeightType computeFarnessDecrement(node v);
    node estimateHighestDecrement();
    node extractQueueTop(std::queue<node> &q);
};
} // namespace GroupClosenessGrowShrinkDetails
} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GROUP_CLOSENESS_GROW_SHRINK_IMPL_HPP_
