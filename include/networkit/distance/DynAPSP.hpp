// no-networkit-format
/*
 * DynAPSP.h
 *
 *  Created on: 12.08.2015
 *      Author: Arie Slobbe, Elisabetta Bergamini
 */

#ifndef NETWORKIT_DISTANCE_DYN_APSP_HPP_
#define NETWORKIT_DISTANCE_DYN_APSP_HPP_

#include <networkit/distance/APSP.hpp>
#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/base/DynAlgorithm.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Dynamic APSP.
 */
class DynAPSP : public APSP, public DynAlgorithm {

public:
    /**
     * Creates the object for @a G.
     *
     * @param G The graph.
     */
    DynAPSP(Graph& G);

    /** initialize distances and Pred by repeatedly running the Dijkstra2 algorithm */
    void run() override;

  /**
  * Updates the pairwise distances after an edge insertions on the graph.
  * Notice: it works only with edge insertions.
  *
  * @param e The edge insertions.
  */
  void update(GraphEvent e) override;

  /**
  * Updates the pairwise distances after a batch of edge insertions on the graph.
  * Notice: it works only with edge insertions.
  *
  * @param batch The batch of edge insertions.
  */
  void updateBatch(const std::vector<GraphEvent>& batch) override;


    /** Returns number of visited pairs */
    count visPairs();

    /**
    * Returns a vector containing a shortest path from node u to node v, and an empty path if u and v are not connected.
    *
    */
    std::vector<node> getPath(node u, node v);

private:

  private:
    count visitedPairs = 0;
};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_DYN_APSP_HPP_
