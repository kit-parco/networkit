// no-networkit-format
/*
 * AdamicAdarIndex.cpp
 *
 *  Created on: 25.03.2015
 *      Author: Kolja Esders
 */

#include <cmath>

#include <networkit/linkprediction/AdamicAdarIndex.hpp>
#include <networkit/linkprediction/NeighborhoodUtility.hpp>

namespace NetworKit {

double AdamicAdarIndex::runImpl(node u, node v) {
  std::vector<node> commonNeighbors = NeighborhoodUtility::getCommonNeighbors(*G, u, v);
  double sum = 0;
  for (node w : commonNeighbors) {
    sum += 1.0 / std::log(static_cast<double>(G->degree(w)));
  }
  return sum;
}

} // namespace NetworKit
