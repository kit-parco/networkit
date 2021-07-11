// no-networkit-format
/*
 * NeighborhoodUtility.hpp
 *
 *  Created on: 06.04.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_NEIGHBORHOOD_UTILITY_HPP_
#define NETWORKIT_LINKPREDICTION_NEIGHBORHOOD_UTILITY_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Provides basic operations on neighborhoods in a given graph.
 */
class NeighborhoodUtility final {
  /**
   * Sorts and returns the neighborhoods of the nodes @a u and @a v in the given graph @a G.
   * @param G Graph to obtain neighborhoods from
   * @param u First node in the graph
   * @param v Second node in the graph
   * @return the sorted nodes regarding the neighborhoods of the nodes @a u and @a v in the given graph @a G.
   */
  static std::pair<std::vector<node>, std::vector<node>> getSortedNeighborhoods(const Graph& G, node u, node v);

public:
  /**
   * Returns the union of the neighboorhoods of @a u and @a v.
   * @param G Graph to obtain neighbors-union from
   * @param u First node
   * @param v Second node
   * @return a vector containing all the nodes in the neighborhood-union of @a u and @a v
   */
  static std::vector<node> getNeighborsUnion(const Graph& G, node u, node v);

  /**
   * Returns a vector containing the node-ids of all common neighbors of @a u and @a v.
   * @param G Graph to obtain common neighbors from
   * @param u First node
   * @param v Second node
   * @return a vector containing the node-ids of all common neighbors of @a u and @a v
   */
  static std::vector<node> getCommonNeighbors(const Graph& G, node u, node v);

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_NEIGHBORHOOD_UTILITY_HPP_
