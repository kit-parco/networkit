// no-networkit-format
/*
 * ROCMetric.hpp
 *
 *  Created on: 14.03.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_ROC_METRIC_HPP_
#define NETWORKIT_LINKPREDICTION_ROC_METRIC_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/linkprediction/EvaluationMetric.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Provides points that define the Receiver Operating Characteristic curve for a given set of predictions.
 * Based on the generated points the area under the curve can be calculated with the trapzoidal rule.
 */
class ROCMetric final: public EvaluationMetric {
  /**
   * Generate the points of the Receiver Operating Characteristic curve regarding the previously set predictions.
   * Note that in the case of multiple y-values mapping to the same x-value the highest (=latest) y-value gets picked.
   * @return a pair of vectors where the first vector contains the false positive rates and the second vector the
   * corresponding true positive rates
   */
  std::pair<std::vector<double>, std::vector<double>> generatePoints() override;

public:
  using EvaluationMetric::EvaluationMetric;

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_ROC_METRIC_HPP_
