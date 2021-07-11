// no-networkit-format
/*
 * EdmondsKarp.hpp
 *
 *  Created on: 11.06.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu), Michael Hamann <michael.hamann@kit.edu>
 */

#ifndef NETWORKIT_FLOW_EDMONDS_KARP_HPP_
#define NETWORKIT_FLOW_EDMONDS_KARP_HPP_

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup flow
 * The EdmondsKarp class implements the maximum flow algorithm by Edmonds and Karp.
 */
class EdmondsKarp final : public Algorithm {
    const Graph *graph;

    node source, sink;

    std::vector<edgeweight> flow;
    edgeweight flowValue;

    /**
     * Performs a breadth-first search on the graph from the source node to find an augmenting path to the sink node respecting the flow values
     * @param residFlow The residual flow in the network.
     * @param pred Used to store the path from the source to the sink.
     * @return The gain in terms of flow.
     */
    edgeweight BFS(std::vector< edgeweight > &residFlow, std::vector< node > &pred) const;

public:
    /**
     * Constructs an instance of the EdmondsKarp algorithm for the given graph, source and sink
     * @param graph The graph.
     * @param source The source node.
     * @param sink The sink node.
     */
     EdmondsKarp(const Graph &graph, node source, node sink);

    /**
     * Computes the maximum flow, executes the EdmondsKarp algorithm.
     */
    void run() override;

    /**
     * Returns the value of the maximum flow from source to sink.
     *
     * @return The maximum flow value
     */
    edgeweight getMaxFlow() const;

    /**
     * Returns the set of the nodes on the source side of the flow/minimum cut.
     *
     * @return The set of nodes that form the (smallest) source side of the flow/minimum cut.
     */
    std::vector<node> getSourceSet() const;

    /**
     * Get the flow value between two nodes @a u and @a v.
     * @warning The running time of this function is linear in the degree of u.
     *
     * @param u The first node
     * @param v The second node
     * @return The flow between node u and v.
     */
    edgeweight getFlow(node u, node v) const;

    /**
     * Get the flow value of an edge.
     *
     * @param eid The id of the edge
     * @return The flow on the edge identified by eid
     */
    edgeweight getFlow(edgeid eid) const {
        assureFinished();
        return flow[eid];
    };

    /**
     * Return a copy of the flow values of all edges.
     * @note Instead of copying all values you can also use the inline function "getFlow(edgeid)" in order to access the values efficiently.
     *
     * @return The flow values of all edges
     */
    const std::vector<edgeweight> &getFlowVector() const;
};

} /* namespace NetworKit */

#endif // NETWORKIT_FLOW_EDMONDS_KARP_HPP_
