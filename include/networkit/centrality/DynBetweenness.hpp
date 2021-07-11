// no-networkit-format
/*
 * DynBetweenness.h
 *
 *  Created on: 12.08.2015
 *      Author: Arie Slobbe, Elisabetta Bergamini
 */

#ifndef NETWORKIT_CENTRALITY_DYN_BETWEENNESS_HPP_
#define NETWORKIT_CENTRALITY_DYN_BETWEENNESS_HPP_

#include <memory>
#include <queue>

#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/auxiliary/PrioQueue.hpp>
#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

class CompareDist
{
public:
    bool operator()(std::pair<double,node> n1,std::pair<double,node> n2) {
        return n1.first > n2.first;
    }
};


/**
 * @ingroup centrality
 * Dynamic APSP.
 */
class DynBetweenness: public Centrality, public DynAlgorithm {

public:
    /**
     * Creates the object for @a G.
     *
     * @param G The graph.
     */
    DynBetweenness(const Graph &G);

  /**
   * Runs static betweenness centrality algorithm on the initial graph.
   */
  void run() override;

  /**
  * Updates the betweenness centralities after an edge insertions on the graph.
  * Notice: it works only with edge insertions.
  *
  * @param e The edge insertions.
  */
  void update(GraphEvent e) override;

  /**
  * Updates the betweenness centralities after a batch of edge insertions on the graph.
  * Notice: it works only with edge insertions.
  *
  * @param batch The batch of edge insertions.
  */
  void updateBatch(const std::vector<GraphEvent>& batch) override;


    /** Returns number of visited pairs */
    count visPairs();

    edgeweight getDistance(node u, node v);
    edgeweight getSigma(node u, node v);

    count numAffectedAPSP();

    count numAffectedDep();

  double getTimeDep();


private:
    void increaseScore(std::vector<bool> & affected, node y, std::priority_queue<std::pair<double, node>, std::vector<std::pair<double,node>>, CompareDist> & Q);
    void decreaseScore(std::vector<bool> & affected, node y, std::priority_queue<std::pair<double, node>, std::vector<std::pair<double,node>>, CompareDist> & Q);
    node u;
    node v;
    edgeweight diameter = 0.0;
    const edgeweight infDist = std::numeric_limits<edgeweight>::max();
    count visitedPairs = 0;
    std::vector<std::vector<edgeweight>> distances;
    std::vector<std::vector<edgeweight>> distancesOld;
    // total number of shortest paths between two nodes
    std::vector<std::vector<edgeweight>> sigma;
    std::vector<std::vector<edgeweight>> sigmaOld;

    count affectedAPSP = 0;
    count affectedDep = 0;
  double timeDep = 0;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_DYN_BETWEENNESS_HPP_
