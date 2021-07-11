// no-networkit-format
/*
 * RandomLinkSampler.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders
 */

#include <networkit/graph/GraphTools.hpp>
#include <networkit/linkprediction/RandomLinkSampler.hpp>

namespace NetworKit {

namespace RandomLinkSampler {

Graph byPercentage(const Graph& G, double trainPercentage) {
  if (trainPercentage < 0 || trainPercentage > 1) {
    throw std::invalid_argument("Given percentage is not in [0, 1].");
  }
  return byCount(G, trainPercentage * G.numberOfEdges());
}

Graph byCount(const Graph& G, count numTrainLinks) {
  if (numTrainLinks > G.numberOfEdges()) {
    throw std::invalid_argument("numTrainLinks > G.numberOfEdges().");
  }
  Graph trainingGraph(G);
  for (count i = 0; i < G.numberOfEdges() - numTrainLinks; ++i) {
    std::pair<node, node> edgeToRemove = GraphTools::randomEdge(trainingGraph);
    trainingGraph.removeEdge(edgeToRemove.first, edgeToRemove.second);
  }
  return trainingGraph;
}

} // namespace RandomLinkSampler

} // namespace NetworKit
