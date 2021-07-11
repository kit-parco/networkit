// no-networkit-format
/*
 * ResourceAllocationIndex.cpp
 *
 *  Created on: 11.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include <networkit/linkprediction/ResourceAllocationIndex.hpp>
#include <networkit/linkprediction/NeighborhoodUtility.hpp>

namespace NetworKit {

double ResourceAllocationIndex::runImpl(node u, node v) {
  std::vector<node> commonNeighbors = NeighborhoodUtility::getCommonNeighbors(*G, u, v);
  double sum = 0;
  for (node w : commonNeighbors) {
    sum += 1.0 / G->degree(w);
  }
  return sum;
}

} // namespace NetworKit
