// no-networkit-format
/*
 * PredictionsSorter.hpp
 *
 *  Created on: 26.04.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_PREDICTIONS_SORTER_HPP_
#define NETWORKIT_LINKPREDICTION_PREDICTIONS_SORTER_HPP_

#include <networkit/linkprediction/LinkPredictor.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Allows the sorting of predictions by score or node-pair.
 */
class PredictionsSorter final {
  /**
   * Comparator used to sort predictions descendingly by score and on equality
   * ascendingly by node-pairs which means e.g. (0, 1) < (1, 1) and (0, 0) < (0, 1).
   */
  struct ScoreComp {
    bool operator()(const LinkPredictor::prediction& a, const LinkPredictor::prediction& b) const {
      return (a.second > b.second) || (a.second == b.second && a.first < b.first);
    }
  };

  static ScoreComp ConcreteScoreComp; //!< Comparator instance for score-based comparison

  /**
   * Comparator used to sort predictions ascendingly by node-pairs.
   */
  struct NodePairComp {
    bool operator()(const LinkPredictor::prediction& a, const LinkPredictor::prediction& b) const {
      return a.first < b.first;
    }
  };

  static NodePairComp ConcreteNodePairComp; //!< Comparator instance for node-pair-based comparison

public:
  /**
   * Sorts the @a predictions descendingly by score.
   * In case there is a tie the node-pairs are used as a tie-breaker by sorting them
   * ascendingly on the first node and on a tie ascendingly by the second node.
   * @param predictions The predictions to sort
   */
  static void sortByScore(std::vector<LinkPredictor::prediction>& predictions);

  /**
   * Sorts the @a predictions ascendingly by node-pair.
   * This means for example (0, 0) < (0, 1) and (1, 1) < (1, 0).
   * @param predictions The predictions to sort
   */
  static void sortByNodePair(std::vector<LinkPredictor::prediction>& predictions);

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_PREDICTIONS_SORTER_HPP_
