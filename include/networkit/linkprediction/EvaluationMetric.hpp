// no-networkit-format
/*
 * EvaluationMetric.hpp
 *
 *  Created on: 17.03.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_EVALUATION_METRIC_HPP_
#define NETWORKIT_LINKPREDICTION_EVALUATION_METRIC_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/linkprediction/LinkPredictor.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Abstract base class for evaluation curves.
 * The evaluation curves are generated based on the predictions calculated
 * by the link predictor and a testGraph to compare against.
 */
class EvaluationMetric {
  /**
   * Generates the points on the evaluation curve. Each metric will implement
   * it's own way to generate those points.
   * Doesn't have to check if testGraph or predictions are existent.
   * @return a pair of X- and Y-vectors (in that order) that can be used for plotting
   */
  virtual std::pair<std::vector<double>, std::vector<double>> generatePoints() = 0;

  /**
   * Helper function that bundles the calculation of all statistical measures including
   * number of TP, FP, TN, FN and also absolute number of positives and negatives.
   */
  void calculateStatisticalMeasures();

  // The following three helper methods generate statistical measures for every threshold
  // based on the given predictions. The names should be self-explanatory.
  void setPositivesAndNegatives();

  void setTrueAndFalsePositives();

  void setTrueAndFalseNegatives();

  /**
   * Sorts the given pair of vectors ascendingly with respect to the values of the first vector.
   * This effectively leads to an sorted pair of vectors where each x-value still maps to the
   * same y-value but the x-values are sorted ascendingly from 0 to 1.
   * @param curve Points to sort
   */
  void sortPointsOfCurve(std::pair<std::vector<double>, std::vector<double>>& curve) const;

protected:
  std::pair<std::vector<double>, std::vector<double>> generatedPoints; //!< Points describing the generated curve. Will be set after a call to getCurve

  const Graph* testGraph; //!< Used to evaluate the binary predictions at the thresholds

  std::vector<LinkPredictor::prediction> predictions; //!< Predictions that should be evaluated

  std::vector<index> thresholds; //!< Indices for the thresholds to use. All node-pairs with an index < thresholds[i] will be regarded as links

  count numPositives; //!< Absolute number of positive instances in the prediction-set

  count numNegatives; //!< Absolute number of negative instances in the prediction-set

  std::vector<count> truePositives; //!< True positives regarding the corresponding threshold

  std::vector<count> falsePositives; //!< False positives regarding the corresponding threshold

  std::vector<count> trueNegatives; //!< True negatives regarding the corresponding threshold

  std::vector<count> falseNegatives; //!< False negatives regarding the corresponding threshold

public:
  EvaluationMetric();

  /**
   *
   * @param testGraph Graph containing the links to use for evaluation
   * @param predictions Dyad-score-pairs whose prediction quality will be evaluated
   */
  explicit EvaluationMetric(const Graph& testGraph);

  /**
   * Default destructor.
   */
  virtual ~EvaluationMetric() = default;

  /**
   * Sets a new graph to use as ground truth for evaluation.
   * Note that this won't reset the most recently calculated curve and as a consequence
   * getAreaUnderCurve() const will still behave as expected by returning the AUC of the most recent curve.
   * @param newTestGraph New graph to use as ground truth
   */
  void setTestGraph(const Graph& newTestGraph);

  /**
   * Returns a pair of X- and Y-vectors describing the evaluation curve generated from the @a predictions.
   * The latest y-value will be used as a tie-breaker in case there are multiple y-values for one x-value.
   * Note that the given number of thresholds (@a numThresholds) is an upper bound for the number of
   * points returned. This is due to the fact that multiple y-values can map to one x-value in which case
   * the tie-breaking behavior described above will intervene.
   * @param predictions Predictions to evaluate
   * @param numThresholds The number of thresholds to use the metric on
   * @return a pair of vectors where the first vectors contains all x-values and the second one contains the corresponding
   * y-value
   */
  virtual std::pair<std::vector<double>, std::vector<double>> getCurve(std::vector<LinkPredictor::prediction> predictions, count numThresholds = 1000);

  /**
   * Returns the area under the given @a curve by using the trapezoidal rule.
   * @param curve Curve whose AUC to determine
   * @return the area under the given curve
   */
  virtual double getAreaUnderCurve(std::pair<std::vector<double>, std::vector<double>> curve) const;

  /**
   * Returns the area under the curve that was most recently calculated by this instance.
   * This implies that getCurve() has to get called beforehand.
   * @return area under the most recently calculated curve
   */
  virtual double getAreaUnderCurve() const;

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_EVALUATION_METRIC_HPP_
