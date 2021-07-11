#ifndef NETWORKIT_CENTRALITY_GROUP_CLOSENESS_LOCAL_SEARCH_HPP_
#define NETWORKIT_CENTRALITY_GROUP_CLOSENESS_LOCAL_SEARCH_HPP_

#include <limits>
#include <memory>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class GroupClosenessLocalSearch final : public Algorithm {
public:
    /**
     * Local search approximation algorithm for Group Closeness Maximization presented in
     * "Group-Harmonic and Group-Closeness Maximization – Approximation and Engineering", Angriman
     * et al., ALENEX 2021. The algorithm evaluates all possible swaps between a vertex in the group
     * and the vertices outside, and performs a swap iff the decrement in farness is at least $$(1 -
     * 1 / (k \\cdot (n - k)))$$, where $$k$$ is the number of vertices in the group. Thus,
     * even in a best-case scenario the time complexity of this algorithm is $$O(n \\cdot k)$$. To
     * keep the number of swaps low, it is recommended to use this algorithm as a refinement step of
     * an already good solution computed by a faster algorithm e.g., greedy (GroupCloseness), or
     * GrowShrink (GroupClosenessGrowShrink). In undirected graphs the approximation ratio is 1/5,
     * on directed graphs it has not been demonstrated.
     *
     * @param G A graph.
     * @param first, last A range that contains the initial group of nodes.
     * @param runGrowShrink Whether or not to run the Grow-Shrink algorithm on the initial group.
     * @param maxIterations Maximum number of swaps allowed. Prevents the algorithm from performing
     * too many swaps by giving up the approximation guarantee.
     */
    template <class InputIt>
    GroupClosenessLocalSearch(const Graph &G, InputIt first, InputIt last,
                              bool runGrowShrink = true,
                              count maxIterations = std::numeric_limits<count>::max())
        : GroupClosenessLocalSearch(G, std::vector<node>(first, last), runGrowShrink,
                                    maxIterations) {}

    GroupClosenessLocalSearch(const Graph &G, const std::vector<node> &group,
                              bool runGrowShrink = true,
                              count maxIterations = std::numeric_limits<count>::max());

    /**
     * Runs the algorithm.
     */
    void run() override;

    /**
     * Returns the computed group.
     */
    std::vector<node> groupMaxCloseness() const;

    /**
     * Returns the number of iterations performed by the algorithm.
     */
    count numberOfIterations() const;

    class GroupClosenessLocalSearchInterface : public Algorithm {
    public:
        template <class InputIt>
        GroupClosenessLocalSearchInterface(InputIt first, InputIt last) : group(first, last) {}

        std::unordered_set<node> group;
        count nIterations;
    };

private:
    const bool weighted;
    // Is always one between GroupClosenessLocalSearchImpl<count/edgeweight>,
    // see implementation.
    std::unique_ptr<GroupClosenessLocalSearchInterface> impl;
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GROUP_CLOSENESS_LOCAL_SEARCH_HPP_
