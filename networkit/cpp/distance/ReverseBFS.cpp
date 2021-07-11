// no-networkit-format
/*
 * BFS.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#include <networkit/distance/ReverseBFS.hpp>
#include <queue>

namespace NetworKit {

ReverseBFS::ReverseBFS(const Graph &G, node source, bool storePaths,
                       bool storeNodesByDistance, node target)
    : SSSP(G, source, storePaths, storeNodesByDistance, target) {}

void ReverseBFS::run() {
  edgeweight infDist = std::numeric_limits<edgeweight>::max();
  count z = G->upperNodeIdBound();
  distances.clear();
  distances.resize(z, infDist);

  std::vector<bool> visited;
  visited.resize(z, false);

  if (storePaths) {
    previous.clear();
    previous.resize(z);
    npaths.clear();
    npaths.resize(z, 0);
    npaths[source] = 1;
  }

  if (storeNodesSortedByDistance) {
    std::vector<node> empty;
    std::swap(nodesSortedByDistance, empty);
  }

  std::queue<node> q;
  q.push(source);
  visited[source] = true;
  distances[source] = 0;
  bool breakWhenFound = (target != none);
  while (!q.empty()) {
    node u = q.front();
    q.pop();

    if (storeNodesSortedByDistance) {
      nodesSortedByDistance.push_back(u);
    }
    if (breakWhenFound && target == u) {
      break;
    }

    // insert untouched neighbors into queue
    G->forInNeighborsOf(u, [&](node v) {
      if (!visited[v]) {
        q.push(v);
        visited[v] = true;
        distances[v] = distances[u] + 1;
        if (storePaths) {
          previous[v] = {u};
          npaths[v] = npaths[u];
        }
      } else if (storePaths && (distances[v] == distances[u] + 1)) {
        previous[v].push_back(u); // additional predecessor
        npaths[v] += npaths[u]; // all the shortest paths to u are also shortest
                                // paths to v now
      }
    });
  }
}

} /* namespace NetworKit */
