// no-networkit-format
/*
 * ROCMetric.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders
 */

#include <networkit/linkprediction/ROCMetric.hpp>

namespace NetworKit {

std::pair<std::vector<double>, std::vector<double>> ROCMetric::generatePoints() {
  if (numPositives == 0) {
    throw std::logic_error("ROC metric is not defined for #positives == 0.");
  } else if (numNegatives == 0) {
    throw std::logic_error("ROC metric is not defined for #negatives == 0.");
  }
  std::pair<std::vector<double>, std::vector<double>> points;
  points.first.reserve(thresholds.size());
  points.second.reserve(thresholds.size());
  for (index i = 0; i < thresholds.size(); ++i) {
    double falsePositiveRatio = 1.0 * falsePositives.at(i) / numNegatives;
    if (!points.first.empty() && points.first.back() == falsePositiveRatio) {
      points.second.pop_back();
    } else {
      points.first.push_back(falsePositiveRatio);
    }
    points.second.push_back(1.0 * truePositives.at(i) / numPositives);
  }
  return points;
}

} // namespace NetworKit
