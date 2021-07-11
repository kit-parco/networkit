// no-networkit-format
/*
 * LinkThresholder.cpp
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders
 */

#include <list>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/linkprediction/LinkThresholder.hpp>
#include <networkit/linkprediction/PredictionsSorter.hpp>

namespace NetworKit {

namespace LinkThresholder {

std::vector<std::pair<node, node>>
byScore(const std::vector<LinkPredictor::prediction> &predictions, double minScore) {
  std::list<LinkPredictor::prediction> predictionsList;
  std::copy(predictions.begin(), predictions.end(), std::back_inserter(predictionsList));
  predictionsList.erase(std::remove_if(predictionsList.begin(), predictionsList.end(),
      [&](const std::pair<std::pair<node, node>, double>& p) { return p.second < minScore; }), predictionsList.end());
  std::vector<std::pair<node, node>> selectedLinks;
  selectedLinks.reserve(predictionsList.size());
  std::transform(predictionsList.begin(), predictionsList.end(), std::back_inserter(selectedLinks),
      [](const LinkPredictor::prediction& p){ return p.first; });
  Aux::Parallel::sort(selectedLinks.begin(), selectedLinks.end());
  return selectedLinks;
}

std::vector<std::pair<node, node>>
byCount(const std::vector<LinkPredictor::prediction> &predictions, count numLinks) {
  if (numLinks > predictions.size()) {
    throw std::invalid_argument("numLinks > predictions.size().");
  }
  auto sortedPredictions = predictions;
  PredictionsSorter::sortByScore(sortedPredictions);
  std::vector<LinkPredictor::prediction> selectedPredictions(sortedPredictions.begin(),
                                                             sortedPredictions.begin() + numLinks);
  std::vector<std::pair<node, node>> selectedLinks;
  selectedLinks.reserve(selectedPredictions.size());
  std::transform(selectedPredictions.begin(), selectedPredictions.end(), std::back_inserter(selectedLinks),
      [](const LinkPredictor::prediction& p){ return p.first; });
  Aux::Parallel::sort(selectedLinks.begin(), selectedLinks.end());
  return selectedLinks;
}

std::vector<std::pair<node, node>>
byPercentage(const std::vector<LinkPredictor::prediction> &predictions, double percentageLinks) {
    if (percentageLinks < 0 || percentageLinks > 1) {
        throw std::invalid_argument("Given percentage is not in [0, 1].");
    }
    return byCount(predictions, percentageLinks * predictions.size());
}

} // namespace LinkThresholder

} // namespace NetworKit
