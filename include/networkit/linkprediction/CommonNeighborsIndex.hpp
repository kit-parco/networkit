// no-networkit-format
/*
 * CommonNeighborsIndex.hpp
 *
 *  Created on: 06.12.2014
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_COMMON_NEIGHBORS_INDEX_HPP_
#define NETWORKIT_LINKPREDICTION_COMMON_NEIGHBORS_INDEX_HPP_

#include <networkit/linkprediction/LinkPredictor.hpp>
#include <networkit/linkprediction/NeighborhoodUtility.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * The CommonNeighborsIndex calculates the number of common
 * neighbors of a node-pair in a given graph.
 */
class CommonNeighborsIndex final : public LinkPredictor {
  /**
   * Returns the number of common neighbors of the given nodes @a u and @a v.
   * @param u First node
   * @param v Second node
   * @return the number of common neighbors of @a u and @a v
   */
  double runImpl(node u, node v) override {
      return NeighborhoodUtility::getCommonNeighbors(*G, u, v).size();
  }

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_COMMON_NEIGHBORS_INDEX_HPP_
