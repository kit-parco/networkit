// no-networkit-format
/*
 * MissingLinksFinder.cpp
 *
 *  Created on: 20.03.2015
 *      Author: Kolja Esders
 */

#include <algorithm>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/linkprediction/MissingLinksFinder.hpp>

namespace NetworKit {

MissingLinksFinder::MissingLinksFinder(const Graph& G) : G(&G) {}

std::vector<std::pair<node, node>> MissingLinksFinder::findAtDistance(count k) {
  std::vector<std::pair<node, node>> missingLinks;
  std::vector<node> nodes;
  nodes.reserve(G->numberOfNodes());
  G->forNodes([&](node u) { nodes.push_back(u); });
  #pragma omp parallel
  {
    std::vector<std::pair<node, node>> missingLinksPrivate;
    #pragma omp for nowait
    for (omp_index i = 0; i < static_cast<omp_index>(nodes.size()); ++i) {
      std::vector<std::pair<node, node>> missingAtU = findFromNode(nodes[i], k);
      // Discard all node-pairs of the form u > v. This removes all duplicates that result from undirected edges.
      missingAtU.erase(std::remove_if(std::begin(missingAtU), std::end(missingAtU),
          [&](std::pair<node, node> p) { return p.first >= p.second; }), std::end(missingAtU));
      missingLinksPrivate.insert(missingLinksPrivate.end(), missingAtU.begin(), missingAtU.end());
    }
    #pragma omp critical
    missingLinks.insert(missingLinks.end(), missingLinksPrivate.begin(), missingLinksPrivate.end());
  }
  DEBUG("Found ", missingLinks.size(), " missing links with distance ", k, ".");
  Aux::Parallel::sort(missingLinks.begin(), missingLinks.end());
  return missingLinks;
}

std::vector<std::pair<node, node>> MissingLinksFinder::findFromNode(node u, count k) {
  std::vector<std::pair<node, node>> missingLinks;
  std::vector<bool> visited;
  visited.resize(G->upperNodeIdBound(), false);
  std::queue<node> q;
  q.push(u);
  visited[u] = true;
  for (count i = 1; i <= k; ++i) {
    std::queue<node> newFound;
    while (!q.empty()) {
      node u = q.front();
      q.pop();
      G->forNeighborsOf(u, [&](node v) {
        if (!visited[v]) {
          newFound.push(v);
          visited[v] = true;
        }
      });
    }
    q = newFound;
  }
  while (!q.empty()) {
      missingLinks.emplace_back(u, q.front());
      q.pop();
  }
  return missingLinks;
}

} // namespace NetworKit
