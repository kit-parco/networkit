/*
 * ApproxSpanningEdge.hpp
 *
 *  Created on: 29.09.2019
 *     Authors: Eugenio Angriman <angrimae@hu-berlin.de>
 *
 */

#ifndef NETWORKIT_CENTRALITY_APPROX_SPANNING_EDGE_HPP_
#define NETWORKIT_CENTRALITY_APPROX_SPANNING_EDGE_HPP_

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class ApproxSpanningEdge final : public Algorithm {

    enum class NodeStatus : unsigned char {
        NOT_IN_COMPONENT,
        IN_SPANNING_TREE,
        IN_RANDOM_WALK,
        NOT_VISITED
    };

public:
    /**
     * Computes an epsilon-approximation of the spanning edge centrality of every edge of the input
     * graph with probability (1 - 1/n), based on "Efficient Algorithms for Spanning Tree
     * Centrality", Hayashi et al., IJCAI, 2016. This implementation also support multi-threading.
     *
     * @param G An undirected graph.
     * @param eps Maximum additive error.
     */
    ApproxSpanningEdge(const Graph &G, double eps = 0.1);

    ~ApproxSpanningEdge() override = default;

    /**
     * Executes the algorithm.
     */
    void run() override;

    /**
     * Return the spanning edge approximation for each edge of the graph.
     *
     * @return Spanning edge approximation for each edge of the input graph.
     */
    std::vector<double> scores() const;

private:
    const Graph &G;
    double eps;
    double delta;
    count nSamples;

    // For each thread, marks which nodes have been visited by the random walk.
    std::vector<std::vector<NodeStatus>> visitedNodes;

    // For each thread, counts how many times each edge appears in a random
    // spanning tree.
    std::vector<std::vector<count>> edgeScores;

    // Sequence of biconnected components.
    std::vector<std::vector<node>> sequences;

    // Parent pointers.
    std::vector<std::vector<node>> parents;

    // Ids of edge connecting the nodes to their parents.
    std::vector<std::vector<edgeid>> parentsEdgeIds;

    // Samples a spanning tree uniformly at random.
    void sampleUST();
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_APPROX_SPANNING_EDGE_HPP_
