/*
 *  GroupClosenessGrowShrink.hpp
 *
 *  Created on: 19.12.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_CENTRALITY_GROUP_CLOSENESS_GROW_SHRINK_HPP_
#define NETWORKIT_CENTRALITY_GROUP_CLOSENESS_GROW_SHRINK_HPP_

#include <memory>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

// pImpl
namespace GroupClosenessGrowShrinkDetails {
template <class>
class GroupClosenessGrowShrinkImpl;
} // namespace GroupClosenessGrowShrinkDetails

class GroupClosenessGrowShrink final : public Algorithm {

public:
    GroupClosenessGrowShrink(const Graph &graph, const std::vector<node> &group,
                             bool extended = false, count insertions = 0,
                             count maxIterations = 100);
    /**
     * Finds a group of nodes with high group closeness centrality. This is the Grow-Shrink
     * algorithm presented in Angriman et al. "Local Search for Group Closeness Maximization on Big
     * Graphs" IEEE BigData 2019. The algorithm takes as input a graph and an arbitrary group of
     * nodes, and improves the group closeness of the given group by performing vertex exchanges.
     *
     * @param G A connected undirected graph.
     * @param first Iterator for first node of initial group of nodes.
     * @param last Iterator for last node of initial group of nodes.
     * @param extended Set this parameter to true for the Extended Grow-Shrink algorithm (i.e.,
     * vertex exchanges are not restricted to only neighbors of the group).
     * @param insertions Number of consecutive node insertions and removal per iteration. Let this
     * parameter to zero to use Diameter(G)/sqrt(k) nodes (where k is the size of the group).
     * @param maxIterations Maximum number of iterations allowed.
     */
    template <class Iter>
    GroupClosenessGrowShrink(const Graph &G, Iter first, Iter last, bool extended = false,
                             count insertions = 0, count maxIterations = 100)
        : GroupClosenessGrowShrink(G, std::vector<node>(first, last), extended, insertions,
                                   maxIterations) {}

    ~GroupClosenessGrowShrink() override;

    /**
     * Runs the algorithm.
     */
    void run() override;

    /**
     * Returns the computed group.
     */
    std::vector<node> groupMaxCloseness() const;

    /**
     * Returns the total number of iterations performed by the algorithm.
     */
    count numberOfIterations() const;

private:
    const Graph *G;
    std::unique_ptr<GroupClosenessGrowShrinkDetails::GroupClosenessGrowShrinkImpl<edgeweight>>
        weightedImpl;
    std::unique_ptr<GroupClosenessGrowShrinkDetails::GroupClosenessGrowShrinkImpl<count>>
        unweightedImpl;
};
} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GROUP_CLOSENESS_GROW_SHRINK_HPP_
