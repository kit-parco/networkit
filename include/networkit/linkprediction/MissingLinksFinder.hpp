// no-networkit-format
/*
 * MissingLinksFinder.hpp
 *
 *  Created on: 20.03.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_MISSING_LINKS_FINDER_HPP_
#define NETWORKIT_LINKPREDICTION_MISSING_LINKS_FINDER_HPP_

#include <utility>

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Allows the user to find missing links in the given graph.
 * The absent links to find are narrowed down by providing a distance
 * that the nodes of the missing links should have.
 * For example in case of distance 2 only node-pairs that would close
 * a triangle in the given graph get returned.
 */
class MissingLinksFinder final {
  const Graph* G; //!< The graph to find missing links in

public:
  /**
   *
   * @param G The graph to find missing links in
   */
  explicit MissingLinksFinder(const Graph& G);

  /**
   * Returns all missing links in the graph that have distance @a k.
   * Note that a distance of @a k actually means that there are @a k different links
   * on the path of the two nodes that are connected through that path.
   * @param k Distance of the absent links
   * @return an ascendingly sorted vector of node-pairs where there is a missing link of distance @a k
   * between the two nodes
   */
  std::vector<std::pair<node, node>> findAtDistance(count k);

  /**
   * Returns all missing links in the graph that have distance @a k and are connected to @a u.
   * Note that a distance of @a k actually means that there are @a k different links
   * on the path of the two nodes that are connected through that path.
   * @param u Node to find missing links from
   * @param k Distance of the absent links
   * @return a vector of node-pairs where there is a missing link of distance @a k
   * between the given node @a u and another node in the graph
   */
  std::vector<std::pair<node, node>> findFromNode(node u, count k);

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_MISSING_LINKS_FINDER_HPP_
