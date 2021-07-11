// no-networkit-format
/*
 * NeighborsMeasureIndex.hpp
 *
 *  Created on: 05.04.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_NEIGHBORS_MEASURE_INDEX_HPP_
#define NETWORKIT_LINKPREDICTION_NEIGHBORS_MEASURE_INDEX_HPP_

#include <networkit/linkprediction/LinkPredictor.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Implementation of the Neighbors Measure Index.
 * This index is also known as Friends Measure and simply returns
 * the number of connections between neighbors of the given nodes u and v.
 */
class NeighborsMeasureIndex final : public LinkPredictor {
  /**
   * Returns the number of connections between neighbors of @a u and @a v.
   * @param u First node
   * @param v Second node
   * @return the number of connections between neighbors of @a u and @a v
   */
  double runImpl(node u, node v) override {
    double neighborConnections = 0;
    G->forNeighborsOf(u, [&](node uNeighbor) {
      G->forNeighborsOf(v, [&](node vNeighbor) {
        // Don't count self-loops
        if (uNeighbor == vNeighbor || G->hasEdge(uNeighbor, vNeighbor)) {
          ++neighborConnections;
        }
      });
    });
    return neighborConnections;
  }

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_NEIGHBORS_MEASURE_INDEX_HPP_
