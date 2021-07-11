// no-networkit-format
/*
 * DYNBETNODE.h
 *
 *  Created on: 10.03.2016
 *      Author: Elisabetta Bergamini
 */

#ifndef NETWORKIT_CENTRALITY_DYN_BETWEENNESS_ONE_NODE_HPP_
#define NETWORKIT_CENTRALITY_DYN_BETWEENNESS_ONE_NODE_HPP_

#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/base/DynAlgorithm.hpp>

namespace NetworKit {

/**
 * @ingroup graph
 * Dynamic betweenness of a single node.
 */
class DynBetweennessOneNode: public DynAlgorithm {

public:
    /**
     * Creates the object for @a G.
     *
     * @param G The graph.
     * @paarm x The node for which we want to compute betweenness.
     */
    DynBetweennessOneNode(Graph& G, node x);

    /** initialize distances and Pred by repeatedly running the Dijkstra2 algorithm */
    void run();

    /**
  * Updates the betweenness centrality of x after an edge insertions on the graph.
  * Notice: it works only with edge insertions.
  *
  * @param e The edge insertions.
  */
  void update(GraphEvent event) override;

    /**
    * Updates the betweenness centrality of x after a batch of edge insertions on the graph.
    * Notice: it works only with edge insertions.
    *
    * @param batch The batch of edge insertions.
    */
    void updateBatch(const std::vector<GraphEvent>& batch) override;

    /** Returns the new betweenness score of node x after the insertion of an edge. Distances and scores of the other nodes are not changed by this function. */
    edgeweight computeScore(GraphEvent event);

    edgeweight getDistance(node u, node v);
    edgeweight getSigma(node u, node v);
    edgeweight getSigmax(node u, node v);
    edgeweight getbcx();

private:
    Graph& G;
    node x;
    // betweenness centrality of node x
    edgeweight bcx = 0;
    const edgeweight infDist = std::numeric_limits<edgeweight>::max();
    const edgeweight epsilon = 0.0000000001; //make sure that no legitimate edge weight is below that.

    std::vector<std::vector<edgeweight>> distances;
    std::vector<std::vector<edgeweight>> distancesOld;
    // total number of shortest paths between two nodes
    std::vector<std::vector<edgeweight>> sigma;
    // number of shortest paths between two nodes that go through node x
    std::vector<std::vector<edgeweight>> sigmax;
    std::vector<std::vector<edgeweight>> sigmaxOld;

    std::vector<node> Pred;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_DYN_BETWEENNESS_ONE_NODE_HPP_
