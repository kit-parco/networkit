// no-networkit-format
/*
 * PrecisionRecallMetric.hpp
 *
 *  Created on: 21.03.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_PRECISION_RECALL_METRIC_HPP_
#define NETWORKIT_LINKPREDICTION_PRECISION_RECALL_METRIC_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/linkprediction/EvaluationMetric.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Provides points that define the Precision-Recall curve for a given set of predictions.
 * Based on the generated points the area under the curve can be calculated with the trapzoidal rule.
 */
class PrecisionRecallMetric final : public EvaluationMetric {
  /**
   * Generates the points for the Precision-Recall curve with respect to the given predictions.
   * The curve assigns every recall-value a corresponding precision as the y-value.
   * In case of a tie regarding multiple y-values for a x-value the smallest (= latest) y-value will be used.
   * @return a pair of vectors where the first vector contains all recall-values and the second vector
   * the corresponding precision-values
   */
  std::pair<std::vector<double>, std::vector<double>> generatePoints() override;

public:
  using EvaluationMetric::EvaluationMetric;

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_PRECISION_RECALL_METRIC_HPP_
