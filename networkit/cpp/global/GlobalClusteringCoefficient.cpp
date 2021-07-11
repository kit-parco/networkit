// no-networkit-format
/*
 * GlobalClusteringCoefficient.cpp
 *
 *  Created on: 12.11.2013
 */

#include <networkit/auxiliary/Random.hpp>
#include <networkit/global/GlobalClusteringCoefficient.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

int uniformRandom(int max) {
  static int offset = 0;
  
  int currentMax = 1;
  int currentValue = 0;
  while(currentMax < max) {
    currentValue = currentValue * RAND_MAX + static_cast<int>(Aux::Random::integer());
    currentMax *= RAND_MAX;
  }
  int value = currentValue % max;
  return offset = (value + offset) % max;
}

unsigned int findIndex(const std::vector<int>& w, int v,
                       unsigned int lowerIdx, unsigned int upperIdx) {
  if(upperIdx - lowerIdx <= 1) {
    return lowerIdx;
  }
  int middleIdx = static_cast<int>((upperIdx + lowerIdx) / 2);
  if(v >= w[middleIdx]) {
    return findIndex(w, v, middleIdx, upperIdx);
  } else {
    return findIndex(w, v, lowerIdx, middleIdx);
  }
}

double GlobalClusteringCoefficient::approximate(const Graph& G, int k) {
  const count n = G.numberOfNodes();
  
  std::vector<int> w(n + 1);
  int sum = 0;
  for(node i = 0; i < n; i++) {
    w[i] = sum;
    sum += static_cast<int>((G.degree(i) * (G.degree(i) - 1)) / 2);
  }
  w[n] = sum;

  int l = 0;
  for(int i = 0; i < k; i++) {
    int r2 = uniformRandom(w[n]);
    node r = findIndex(w, r2, 0, n);
    node u = GraphTools::randomNeighbor(G, r);
    node w;
    do {
      w = GraphTools::randomNeighbor(G, r);
    } while (w == u);
    if(G.hasEdge(u, w)) {
      l++;
    }
  }

  return (double)l / (double)k;
}

} /* namespace NetworKit */
