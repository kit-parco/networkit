// no-networkit-format
/*
 * KatzIndex.hpp
 *
 *  Created on: 30.01.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_KATZ_INDEX_HPP_
#define NETWORKIT_LINKPREDICTION_KATZ_INDEX_HPP_

#include <unordered_map>

#include <networkit/linkprediction/LinkPredictor.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Implementation of the Katz index.
 * Katz index assigns a pair of nodes a similarity score
 * that is based on the sum of the weighted number of paths of length l
 * where l is smaller than a given limit.
 */
class KatzIndex final : public LinkPredictor {
  count maxPathLength; //!< Maximal length of the paths to consider

  double dampingValue; //!< Damping factor in (0,1) used to exponentially damp every addend of the sum

  node lastStartNode; //!< Node at which the algorithm started at the last call

  std::unordered_map<node, double> lastScores; //!< Stores the last generated scores from the lastStartNode

  std::vector<double> dampingFactors; //!< Stores precalculated damping factors to increase performance

  // Helper method used to access the score for a given node-pair. Checks which of the given nodes
  // was used as the starting node and uses the other node to access the last scores generated.
  // Defaults to 0.0 if no score could be found.
  double getScore(node u, node v) const;

  /**
   * Returns the similarity score for the given node-pair based on the Katz index specified during construction.
   * The algorithm considers all paths starting at the node with the smaller degree except the algorithm
   * started at the other node at the last call.
   * @param u First node
   * @param v Second node
   * @return the similarity score of the given node-pair calculated by the specified Katz index
   */
  double runImpl(node u, node v) override;

  /**
   * Calculated the damping factors for every path length smaller or equal to maxPathLength.
   * This can be used to cache damping factors to increase performance through reuse.
   * The results will be stored in dampingFactors and dampingFactors[0] is always 1.
   */
  void calcDampingFactors();

public:
  /**
   * @param maxPathLength Maximal length of the paths to consider
   * @param dampingValue Used to exponentially damp every addend of the sum. Should be in (0, 1]
   */
  explicit KatzIndex(count maxPathLength = 5, double dampingValue = 0.005);

  /**
   * @param G The graph to operate on
   * @param maxPathLength Maximal length of the paths to consider
   * @param dampingValue Used to exponentially damp every addend of the sum. Should be in (0, 1]
   */
  explicit KatzIndex(const Graph& G, count maxPathLength = 5, double dampingValue = 0.005);

  // Overriding this method is necessary as the implementation of the Katz index makes use
  // of caching. This makes run() not thread-safe. To still achieve performance gains
  // we split the nodePairs into subsets and create a new Katz instance for every subset.
  std::vector<LinkPredictor::prediction> runOn(std::vector<std::pair<node, node>> nodePairs) override;

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_KATZ_INDEX_HPP_
