/*
 *  GroupClosenessLocalSwaps.hpp
 *
 *  Created on: 19.12.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_CENTRALITY_GROUP_CLOSENESS_LOCAL_SWAPS_HPP_
#define NETWORKIT_CENTRALITY_GROUP_CLOSENESS_LOCAL_SWAPS_HPP_

#ifdef __AVX2__
#include <immintrin.h>
#include <networkit/auxiliary/AlignedAllocator.hpp>
#endif // __AVX2__

#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class GroupClosenessLocalSwaps final : public Algorithm {

    static constexpr float maxInt16 = static_cast<float>(std::numeric_limits<uint16_t>::max());
    // Number of random numbers generated at once
    static constexpr count K = 16;

public:
    GroupClosenessLocalSwaps(const Graph &G, const std::vector<node> &group, count maxSwaps = 100);

    /**
     * Finds a group of nodes with high group closeness centrality. This is the LS-restrict
     * algorithm presented in Angriman et al. "Local Search for Group Closeness Maximization on Big
     * Graphs" IEEE BigData 2019. The algorithm takes as input a graph and an arbitrary group of
     * nodes, and improves the group closeness of the given group by performing vertex exchanges.
     *
     * @param G A connected, undirected, and unweighted graph.
     * @param first, last A range that contains the initial group of nodes.
     * @param maxSwaps Maximum number of vertex exchanges allowed.
     */
    template <class InputIt>
    GroupClosenessLocalSwaps(const Graph &graph, InputIt first, InputIt last, count maxSwaps = 100)
        : G(&graph), group(first, last), maxSwaps(maxSwaps) {

        if (G->isDirected())
            throw std::runtime_error("Error, this algorithm does not support directed graphs.");

        if (group.empty())
            throw std::runtime_error("Error, empty group.");

        if (G->isWeighted())
            WARN("This algorithm does not support edge Weights, they will be ignored.");
    }

    /**
     * Runs the algorithm.
     */
    void run() override;

    /**
     * Returns the computed group.
     */
    std::vector<node> groupMaxCloseness() const;

    /**
     * Returns the total number of vertex exchanges performed by the algorithm.
     */
    count numberOfSwaps() const;

private:
    const Graph *G;
    std::vector<node> group, stack;
    std::vector<uint32_t> distance, sumOfMins;
    std::vector<bool> gamma, visited;
    std::unordered_map<node, index> idxMap;
    std::vector<int64_t> farness, farnessDecrease;

    count maxSwaps, totalSwaps, stackSize;

    void init();
    void bfsFromGroup();
    bool findAndSwap();
    node estimateHighestDecrease();
    void initRandomVector();
    int64_t computeFarnessDecrease(node u);
    void resetGamma(node x, index idx);

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
};
} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GROUP_CLOSENESS_LOCAL_SWAPS_HPP_
