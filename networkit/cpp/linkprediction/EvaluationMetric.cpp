// no-networkit-format
/*
 * EvaluationMetric.cpp
 *
 *  Created on: 17.03.2015
 *      Author: Kolja Esders
 */

#include <algorithm>
#include <numeric>
#include <set>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/linkprediction/EvaluationMetric.hpp>
#include <networkit/linkprediction/PredictionsSorter.hpp>

namespace NetworKit {

EvaluationMetric::EvaluationMetric() : testGraph(nullptr) {}

EvaluationMetric::EvaluationMetric(const Graph& testGraph) : testGraph(&testGraph) {}

void EvaluationMetric::setTestGraph(const Graph& newTestGraph) {
  testGraph = &newTestGraph;
}

std::pair<std::vector<double>, std::vector<double>> EvaluationMetric::getCurve(std::vector<LinkPredictor::prediction> predictions, count numThresholds) {
  if (testGraph == nullptr) {
    throw std::logic_error("Set testGraph first.");
  } else if (predictions.empty()) {
      throw std::logic_error("predictions.size() == 0");
  } else if (numThresholds < 2) {
      throw std::invalid_argument("numThresholds < 2: At least 2 thresholds needed for curve.");
  }
  // Don't overshoot with the number of thresholds being greater than the number of predictions + 1.
  if (predictions.size() + 1 < numThresholds || numThresholds == 0) {
    numThresholds = predictions.size() + 1;
  }
  std::set<index> thresholdSet;
  count numPredictions = predictions.size();
  for (index i = 0; i < numThresholds; ++i) {
    // Percentile calculation through nearest rank method.
    // This calculation is numerically instable. This is why we use a set to make sure that
    // we obtain pairwise different thresholds.
    thresholdSet.insert(std::ceil(numPredictions * (1.0 * i / static_cast<double>(numThresholds - 1))));
  }
  thresholds.assign(thresholdSet.begin(), thresholdSet.end());
  // The extraction of statistical measures requires sorted predictions
  PredictionsSorter::sortByScore(predictions);
  this->predictions = predictions;
  calculateStatisticalMeasures();
  generatedPoints = generatePoints();
  return generatedPoints;
}

double EvaluationMetric::getAreaUnderCurve(std::pair<std::vector<double>, std::vector<double>> curve) const {
  if (curve.first.size() < 2) {
    throw std::invalid_argument("At least 2 points needed for AUC.");
  } else if (curve.first.size() != curve.second.size()) {
    throw std::invalid_argument("X- and Y-vector differ in size().");
  }
  // Sort points so that x-values start at zero.
  sortPointsOfCurve(curve);
  double AUC = 0;
  for (index i = 0; i < curve.first.size() - 1; ++i) {
    // Trapezoid rule
    AUC += 0.5 * (curve.first[i + 1] - curve.first[i]) * (curve.second[i] + curve.second[i + 1]);
  }
  return AUC;
}

double EvaluationMetric::getAreaUnderCurve() const {
    if (generatedPoints.first.empty()) {
        throw std::logic_error("Call getCurve first or use getAreaUnderCurve(curve).");
    }
  return getAreaUnderCurve(generatedPoints);
}

void EvaluationMetric::calculateStatisticalMeasures() {
  #pragma omp parallel sections
  {
    #pragma omp section
    setTrueAndFalsePositives();
    #pragma omp section
    setTrueAndFalseNegatives();
    #pragma omp section
    setPositivesAndNegatives();
  }
}

void EvaluationMetric::setTrueAndFalsePositives() {
  count tmpTruePositives = 0;
  count tmpFalsePositives = 0;
  truePositives.clear();
  falsePositives.clear();
  std::vector<index>::iterator thresholdIt = thresholds.begin();
  for (index i = 0; i < predictions.size(); ++i) {
    // Check in the beginning to catch threshold = 0 as well.
    if (thresholdIt != thresholds.end() && i == *thresholdIt) {
      truePositives.push_back(tmpTruePositives);
      falsePositives.push_back(tmpFalsePositives);
      ++thresholdIt;
    }
    bool hasEdge = testGraph->hasEdge(predictions[i].first.first, predictions[i].first.second);
    if (hasEdge) {
      tmpTruePositives++;
    } else {
      tmpFalsePositives++;
    }
  }
  // If there is a threshold at predictions.size() add it as well
  if (thresholdIt != thresholds.end()) {
    truePositives.push_back(tmpTruePositives);
    falsePositives.push_back(tmpFalsePositives);
  }
}

void EvaluationMetric::setTrueAndFalseNegatives() {
  count tmpTrueNegatives = 0;
  count tmpFalseNegatives = 0;
  trueNegatives.clear();
  falseNegatives.clear();
  std::vector<index>::reverse_iterator thresholdIt = thresholds.rbegin();
  for (index i = predictions.size(); i > 0; --i) {
    // Check in the beginning to catch threshold = predictions.size().
    if (thresholdIt != thresholds.rend() && i == *thresholdIt) {
      trueNegatives.push_back(tmpTrueNegatives);
      falseNegatives.push_back(tmpFalseNegatives);
      if (*thresholdIt != 0)
        ++thresholdIt;
    }
    bool hasEdge = testGraph->hasEdge(predictions[i - 1].first.first, predictions[i - 1].first.second);
    if (hasEdge) {
      tmpFalseNegatives++;
    } else {
      tmpTrueNegatives++;
    }
  }
  // If there is a threshold at 0 add it as well
  if (thresholdIt != thresholds.rend()) {
    trueNegatives.push_back(tmpTrueNegatives);
    falseNegatives.push_back(tmpFalseNegatives);
  }
  // We need to reverse so that TN/FN and TP/FP have the same index for a given threshold.
  // As we have started at the bottom of the predictions we need to reverse our results since
  // the TP/FP calculation started at the top.
  std::reverse(trueNegatives.begin(), trueNegatives.end());
  std::reverse(falseNegatives.begin(), falseNegatives.end());
}

void EvaluationMetric::setPositivesAndNegatives() {
  numPositives = 0;
  #pragma omp parallel for
  for (omp_index i = 0; i < static_cast<omp_index>(predictions.size()); ++i) {
    if (testGraph->hasEdge(predictions[i].first.first, predictions[i].first.second)) {
      #pragma omp atomic
      numPositives++;
    }
  }
  numNegatives = predictions.size() - numPositives;
}

void EvaluationMetric::sortPointsOfCurve(std::pair<std::vector<double>, std::vector<double>>& curve) const {
  // Create a permutation that you would use in sorting the x-values ascendingly
  std::vector<int> permutation(curve.first.size());
  std::iota(permutation.begin(), permutation.end(), 0);
  Aux::Parallel::sort(permutation.begin(), permutation.end(), [&](int i, int j){ return curve.first[i] < curve.first[j]; });
  // Actually sort x-vector
  Aux::Parallel::sort(curve.first.begin(), curve.first.end());
  // Apply the permutation to the y-vector
  std::vector<double> sortedY(permutation.size());
  std::transform(permutation.begin(), permutation.end(), sortedY.begin(), [&](int i){ return curve.second[i]; });
  curve.second = sortedY;
}

} // namespace NetworKit
