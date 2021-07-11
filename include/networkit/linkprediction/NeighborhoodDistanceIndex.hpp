// no-networkit-format
/*
 * NeighborhoodDistanceIndex.hpp
 *
 *  Created on: 24.06.2013
 *      Authors: cls, Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_NEIGHBORHOOD_DISTANCE_INDEX_HPP_
#define NETWORKIT_LINKPREDICTION_NEIGHBORHOOD_DISTANCE_INDEX_HPP_

#include <math.h>

#include <networkit/linkprediction/LinkPredictor.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Assigns a distance value to pairs of nodes according to the
 * overlap of their neighborhoods.
 */
class NeighborhoodDistanceIndex final : public LinkPredictor {
  /**
   * Returns the Neighborhood Distance index for the given node-pair (@a u, @a v).
   * @param u First node
   * @param v Second node
   * @return the Neighborhood Distance index for the given node-pair (@a u, @a v)
   */
  double runImpl(node u, node v) override {
    count uNeighborhood = G->degree(u);
    count vNeighborhood = G->degree(v);
    count intersection = NeighborhoodUtility::getCommonNeighbors(*G, u, v).size();
    return ((double)intersection) / (sqrt(uNeighborhood * vNeighborhood));
  }

public:
  using LinkPredictor::LinkPredictor;

};

} /* namespace NetworKit */
#endif // NETWORKIT_LINKPREDICTION_NEIGHBORHOOD_DISTANCE_INDEX_HPP_
