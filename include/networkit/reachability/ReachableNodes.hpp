
#ifndef NETWORKIT_REACHABILITY_REACHABLE_NODES_HPP_
#define NETWORKIT_REACHABILITY_REACHABLE_NODES_HPP_

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup reachability
 */
class ReachableNodes final : public Algorithm {
public:
    /**
     * Determines or estimates the number of reachable nodes from each node in the graph.
     *
     * @param G The graph.
     * @param exact Whether or not to compute the number of reachable nodes exactly. Only used for
     * directed graphs, on undirected graphs the number of reachable nodes from every node can be
     * computed in linear time.
     */
    ReachableNodes(const Graph &G, bool exact = true);

    /**
     * Runs the algorithm.
     */
    void run() override;

    /**
     * Returns the number of reachable nodes from the given node @a u. Only available if @a exact is
     * true.
     *
     * @param u A node.
     * @returns The number of nodes reachable from @a u.
     */
    count numberOfReachableNodes(node u) const {
        assureFinished();
        if (!exact)
            throw std::runtime_error("The number of nodes is not computed exactly, run the "
                                     "algorithm with exact = true.");
        return reachableLB[u];
    }

    /**
     * Returns a lower bound of the number of reachable nodes from the given node @a u.
     *
     * @param u A node.
     * @returns Lower bound of number of nodes reachable from @a u.
     */
    count numberOfReachableNodesLB(node u) const { return reachableLB[u]; }

    /**
     * Returns an upper bound of the number of reachable nodes from the given node @a u.
     *
     * @param u A node.
     * @returns Upper bound of number of nodes reachable from @a u.
     */
    count numberOfReachableNodesUB(node u) const {
        assureFinished();
        if (G->isDirected())
            return exact ? reachableLB[u] : reachableUB[u];
        else
            return reachableLB[u];
    }

    /*
     * Whether or not to compute the reachable nodes exactly.
     */
    bool exact;

private:
    const Graph *G;
    std::vector<count> reachableLB, reachableUB;

    void runDirected();
    void runUndirected();
};
} // namespace NetworKit

#endif // NETWORKIT_REACHABILITY_REACHABLE_NODES_HPP_
