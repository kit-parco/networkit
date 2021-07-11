// no-networkit-format
/*
 * SameCommunityIndex.hpp
 *
 *  Created on: 07.04.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_SAME_COMMUNITY_INDEX_HPP_
#define NETWORKIT_LINKPREDICTION_SAME_COMMUNITY_INDEX_HPP_

#include <networkit/linkprediction/LinkPredictor.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Index to determine whether two nodes are in the same community.
 */
class SameCommunityIndex final : public LinkPredictor {
  Partition communities; //!< The communities of the current graph

  /**
   * Returns 1 if the given nodes @a u and @a v are in the same community, 0 otherwise.
   * @param u First node
   * @param v Second node
   * @return 1 if the given nodes @a u and @a v are in the same community, 0 otherwise
   */
  double runImpl(node u, node v) override;

public:
  SameCommunityIndex();

  /**
   *
   * @param G The graph to work on
   */
  explicit SameCommunityIndex(const Graph& G);

  /**
   * Sets the graph to work on.
   * @param newGraph The graph to work on
   */
  void setGraph(const Graph& newGraph) override;

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_SAME_COMMUNITY_INDEX_HPP_
