// no-networkit-format
/*
 * AdjustedRandIndex.cpp
 *
 *  Created on: 11.04.2015
 *      Author: Kolja Esders
 */

#include <set>

#include <networkit/linkprediction/AdjustedRandIndex.hpp>

namespace NetworKit {

double AdjustedRandIndex::runImpl(node u, node v) {
  std::set<node> uNeighbors;
  std::set<node> vNeighbors;

  G->forNeighborsOf(u, [&](node z) {
    uNeighbors.insert(z);
  });
  G->forNeighborsOf(v, [&](node z) {
    vNeighbors.insert(z);
  });
  std::vector<node> commonNeighbors;
  std::set_intersection(uNeighbors.begin(), uNeighbors.end(), vNeighbors.begin(), vNeighbors.end(), std::back_inserter(commonNeighbors));
  std::vector<node> unionNeighbors;
  std::set_union(uNeighbors.begin(), uNeighbors.end(), vNeighbors.begin(), vNeighbors.end(), std::back_inserter(unionNeighbors));
  std::vector<node> uDifference;
  std::set_union(uNeighbors.begin(), uNeighbors.end(), vNeighbors.begin(), vNeighbors.end(), std::back_inserter(uDifference));
  std::vector<node> vDifference;
  std::set_union(vNeighbors.begin(), vNeighbors.end(), uNeighbors.begin(), uNeighbors.end(), std::back_inserter(vDifference));

  double a = commonNeighbors.size();
  double b = uDifference.size();
  double c = vDifference.size();
  double d = static_cast<double>(G->numberOfNodes() - unionNeighbors.size());
  double ad = a*d;
  // Make sure to not divide by zero
  double denominator = (a*b + a*c + 2*ad + b*b + b*d + c*c + c*d);
  return denominator == 0 ? 0 : (2*(ad - b*c)) / denominator;
}

} // namespace NetworKit
