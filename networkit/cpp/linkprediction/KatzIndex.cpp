// no-networkit-format
/*
 * KatzIndex.cpp
 *
 *  Created on: 30.01.2015
 *      Author: Kolja Esders
 */

#include <list>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/linkprediction/KatzIndex.hpp>
#include <networkit/linkprediction/PredictionsSorter.hpp>

namespace NetworKit {

KatzIndex::KatzIndex(count maxPathLength, double dampingValue)
    : maxPathLength(maxPathLength), dampingValue(dampingValue) {
  calcDampingFactors();
}

KatzIndex::KatzIndex(const Graph& G, count maxPathLength, double dampingValue)
    : LinkPredictor(G), maxPathLength(maxPathLength), dampingValue(dampingValue) {
  calcDampingFactors();
}

double KatzIndex::getScore(node u, node v) const {
  node endNode = lastStartNode == u ? v : u;
  if (lastScores.find(endNode) == lastScores.end()) {
    return 0;
  }
  return lastScores.at(endNode);
}

double KatzIndex::runImpl(node u, node v) {
  if (validCache && (lastStartNode == u || lastStartNode == v)) {
    return getScore(u, v);
  }
  lastScores.clear();
  validCache = false;
  std::list<node> toProcess;
  // Start at the node with less neighbors to potentially increase performance
  lastStartNode = G->degree(u) > G->degree(v) ? v : u;
  toProcess.push_front(lastStartNode);
  for (index pathLength = 1; pathLength <= maxPathLength; ++pathLength) {
    std::unordered_map<node, count> hits;
    for (std::list<node>::const_iterator it = toProcess.begin(); it != toProcess.end(); ++it) {
      const node current = *it;
      G->forNeighborsOf(current, [&](node neighbor) {
        hits[neighbor] += 1;
      });
    }
    // Add found nodes to the todo-list for the next round and update scores
    toProcess.clear();
    for (auto kv : hits) {
      lastScores[kv.first] += dampingFactors[pathLength] * kv.second;
      toProcess.push_front(kv.first);
    }
  }
  validCache = true;
  return getScore(u, v);
}

std::vector<LinkPredictor::prediction> KatzIndex::runOn(std::vector<std::pair<node, node>> nodePairs) {
  // Make sure the nodePairs are sorted. This will make use of the caching of the Katz index
  // and will exploit locality in the form of cpu caching as well.
  Aux::Parallel::sort(nodePairs.begin(), nodePairs.end());
  std::vector<prediction> predictions(nodePairs.size());
  KatzIndex katz(*G, maxPathLength, dampingValue);
  #pragma omp parallel
  {
    KatzIndex katz(*G, maxPathLength, dampingValue);
    #pragma omp for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(nodePairs.size()); ++i) {
      predictions[i] = std::make_pair(nodePairs[i], katz.run(nodePairs[i].first, nodePairs[i].second));
    }
  }
  PredictionsSorter::sortByNodePair(predictions);
  return predictions;
}

void KatzIndex::calcDampingFactors() {
  dampingFactors.resize(maxPathLength + 1);
  dampingFactors[0] = 1;
  for (count i = 1; i <= maxPathLength; ++i) {
    dampingFactors[i] = std::pow(dampingValue, i);
  }
}

} // namespace NetworKit
