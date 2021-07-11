// no-networkit-format
/*
 * AffectedNodes.cpp
 * Author: paddya
 *
 * Created on 17. Dezember 2016, 16:35
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/distance/AffectedNodes.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/ReverseBFS.hpp>
#include <networkit/distance/SSSP.hpp>

namespace NetworKit {

AffectedNodes::AffectedNodes(const Graph &G, const GraphEvent &event)
    : G(G), event(event),
      distances(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
      improvements(G.upperNodeIdBound(), 0) {}

void AffectedNodes::run() {

  if (event.type == GraphEvent::EDGE_ADDITION) {
    addedEdge();
  } else if (event.type == GraphEvent::EDGE_REMOVAL) {
    removedEdge();
  }

  hasRun = true;
}

void AffectedNodes::addedEdge() {

  node u = event.u;
  node v = event.v;

  std::vector<node> uniqueAffectedNodes;
  if (!G.isDirected()) {
    // Pre-compute the distances from u and v to each other node in the graph
    std::vector<edgeweight> distancesU;
    std::vector<edgeweight> distancesV;

    std::pair<std::vector<node>, std::vector<edgeweight>> affectedU;
    std::pair<std::vector<node>, std::vector<edgeweight>> affectedV;

    std::vector<node> affectedNodesU;
    std::vector<node> affectedNodesV;

#pragma omp parallel sections
    {
#pragma omp section
      {
        distancesU = bfsWithoutStartNeighbor(u, v);
        affectedU = getAffectedNodes(u, distancesU);

        affectedNodesU = affectedU.first;
        std::sort(affectedNodesU.begin(), affectedNodesU.end());
      }
#pragma omp section
      {
        distancesV = bfsWithoutStartNeighbor(v, u);
        affectedV = getAffectedNodes(v, distancesV);
        affectedNodesV = affectedV.first;
        std::sort(affectedNodesV.begin(), affectedNodesV.end());
      }
    }

    uniqueAffectedNodes = affectedNodesU;

    uniqueAffectedNodes.insert(uniqueAffectedNodes.end(),
                               affectedNodesV.begin(), affectedNodesV.end());

    for (node w : uniqueAffectedNodes) {
      distances[w] = std::min(distancesU[w], distancesV[w]);
    }

    std::vector<edgeweight> newDistancesU = affectedU.second;
    std::vector<edgeweight> newDistancesV = affectedV.second;

    // complete distances first
    G.parallelForNodes([&](node v) {
      newDistancesU[v] = std::min(newDistancesU[v], distancesU[v]);
      newDistancesV[v] = std::min(newDistancesV[v], distancesV[v]);
    });

    // Compute distance distributions
    std::vector<count> prevNumNodesOnLevelFromU(G.upperNodeIdBound(), 0);
    std::vector<count> prevNumNodesOnLevelFromV(G.upperNodeIdBound(), 0);
    std::vector<count> numNodesOnLevelFromU(G.upperNodeIdBound(), 0);
    std::vector<count> numNodesOnLevelFromV(G.upperNodeIdBound(), 0);

    edgeweight infDist = std::numeric_limits<edgeweight>::max();

    count prevNumLevelsU = 0;
    count prevNumLevelsV = 0;
    count numLevelsU = 0;
    count numLevelsV = 0;

    G.forNodes([&](node v) {
      if (distancesU[v] < infDist) {
        prevNumNodesOnLevelFromU[distancesU[v]] += 1;

        if (distancesU[v] > prevNumLevelsU) {
          prevNumLevelsU = distancesU[v];
        }
      }

      if (distancesV[v] < infDist) {
        prevNumNodesOnLevelFromV[distancesV[v]] += 1;

        if (distancesV[v] > prevNumLevelsV) {
          prevNumLevelsV = distancesV[v];
        }
      }

      if (newDistancesU[v] < infDist) {
        numNodesOnLevelFromU[newDistancesU[v]] += 1;

        if (newDistancesU[v] > numLevelsU) {
          numLevelsU = newDistancesU[v];
        }
      }

      if (newDistancesV[v] < infDist) {
        numNodesOnLevelFromV[newDistancesV[v]] += 1;

        if (newDistancesV[v] > numLevelsV) {
          numLevelsV = newDistancesV[v];
        }
      }
    });

    DEBUG("Num levels: ", prevNumLevelsU, " / ", prevNumLevelsV, " / ",
          numLevelsU, " / ", numLevelsV);

    closenessU = 0;
    for (count i = 1; i <= numLevelsU; i++) {
      closenessU += numNodesOnLevelFromU[i] * (1.0 / i);
    }

    closenessV = 0;
    for (count i = 1; i <= numLevelsV; i++) {
      closenessV += numNodesOnLevelFromV[i] * (1.0 / i);
    }

    // For each level from u, compute the maximum improvement
    std::vector<edgeweight> levelImprovementU(numLevelsU + 1, 0);
    std::vector<edgeweight> levelImprovementV(numLevelsV + 1, 0);

    for (index i = 1.0; i <= numLevelsV; i++) {
      for (index j = 1.0; j <= numLevelsU; j++) {
          levelImprovementV[i] +=
              (numNodesOnLevelFromU[j] * 1.0 / static_cast<double>(i + j))
              - (prevNumNodesOnLevelFromU[j] * 1.0 / static_cast<double>(i + j));
      }
    }

    for (index i = 1.0; i <= numLevelsU; i++) {
      for (index j = 1.0; j <= numLevelsV; j++) {
          levelImprovementU[i] +=
              (numNodesOnLevelFromV[j] * 1.0 / static_cast<double>(i + j))
              - (prevNumNodesOnLevelFromV[j] * 1.0 / static_cast<double>(i + j));
      }
    }

    // Update improvements for each affected node
    for (node v : affectedNodesU) {
      improvements[v] = levelImprovementU[newDistancesV[v]];
    }

    for (node v : affectedNodesV) {
      improvements[v] = levelImprovementV[newDistancesU[v]];
    }

  } else {
    // If the graph is directed, we only need two reverse searches
    // from v

    std::vector<edgeweight> distancesU =
        reverseBfsWithoutStartNeighbor(u, none);
    std::vector<edgeweight> distancesV = reverseBfsWithoutStartNeighbor(v, u);

    std::pair<std::vector<node>, std::vector<edgeweight>> affectedV =
        getAffectedNodesBackwards(v, distancesV);

    uniqueAffectedNodes = affectedV.first;

    std::vector<edgeweight> newDistancesV = affectedV.second;

    for (node w : uniqueAffectedNodes) {
      distances[w] = std::min(distancesV[w], distancesU[w]);
    }

    std::vector<edgeweight> oldDistancesU = bfsWithoutStartNeighbor(u, v);

    std::vector<edgeweight> newDistancesU =
        getAffectedNodes(u, oldDistancesU).second;

    // complete distances first
    G.parallelForNodes([&](node v) {
      newDistancesU[v] = std::min(newDistancesU[v], oldDistancesU[v]);
    });

    // Compute distance distributions
    std::vector<count> prevNumNodesOnLevelFromU(G.upperNodeIdBound(), 0);
    std::vector<count> numNodesOnLevelFromU(G.upperNodeIdBound(), 0);

    edgeweight infDist = std::numeric_limits<edgeweight>::max();

    count numLevelsU = 0;
    count numLevelsV = 0;

    G.forNodes([&](node v) {
      if (oldDistancesU[v] < infDist) {
        prevNumNodesOnLevelFromU[oldDistancesU[v]] += 1;
      }

      if (newDistancesU[v] < infDist) {
        numNodesOnLevelFromU[newDistancesU[v]] += 1;

        if (newDistancesU[v] > numLevelsU) {
          numLevelsU = newDistancesU[v];
        }
      }

      if (newDistancesV[v] < infDist) {
        if (newDistancesV[v] > numLevelsV) {
          numLevelsV = newDistancesV[v];
        }
      }
    });

    closenessU = 0;
    for (count i = 1; i <= numLevelsU; i++) {
      closenessU += numNodesOnLevelFromU[i] * (1.0 / i);
    }

    // For each level from u, compute the maximum improvement
    std::vector<edgeweight> levelImprovementV(numLevelsV + 1, 0);

    for (index i = 1.0; i <= numLevelsV; i++) {
      for (index j = 1.0; j <= numLevelsU; j++) {
          levelImprovementV[i] +=
              (numNodesOnLevelFromU[j] * 1.0 / static_cast<double>(i + j))
              - (prevNumNodesOnLevelFromU[j] * 1.0 / static_cast<double>(i + j));
      }
    }

    for (node v : uniqueAffectedNodes) {
      improvements[v] = levelImprovementV[newDistancesV[v] - 1];
    }
  }

  nodes = std::move(uniqueAffectedNodes);
}

void AffectedNodes::removedEdge() {

  node u = event.u;
  node v = event.v;

  std::vector<node> uniqueAffectedNodes;
  if (!G.isDirected()) {

    BFS bfsU(G, u);
    bfsU.run();
    auto distancesU = bfsU.getDistances();

    BFS bfsV(G, v);
    bfsV.run();
    auto distancesV = bfsV.getDistances();

    // Run a second BFS from u and v but abort if the distance to the node has
    // not decreased
    std::pair<std::vector<node>, std::vector<edgeweight>> affectedU =
        getAffectedNodes(u, distancesU, v);
    std::pair<std::vector<node>, std::vector<edgeweight>> affectedV =
        getAffectedNodes(v, distancesV, u);

    std::vector<node> affectedNodesU = affectedU.first;
    std::vector<node> affectedNodesV = affectedV.first;

    std::sort(affectedNodesU.begin(), affectedNodesU.end());
    std::sort(affectedNodesV.begin(), affectedNodesV.end());

    affectedNodesU.insert(affectedNodesU.end(), affectedNodesV.begin(),
                          affectedNodesV.end());
    uniqueAffectedNodes = std::move(affectedNodesU);

    for (node w : uniqueAffectedNodes) {
      distances[w] = std::min(distancesU[w], distancesV[w]);
    }

  } else {
    // If the graph is directed, we only need two reverse searches
    // from v

    ReverseBFS bfsV(G, v);
    bfsV.run();
    auto distancesV = bfsV.getDistances();

    // G.addEdge(u, v);

    uniqueAffectedNodes = getAffectedNodesBackwards(v, distancesV, u).first;

    for (node w : uniqueAffectedNodes) {
      distances[w] = distancesV[w];
    }
  }

  nodes = std::move(uniqueAffectedNodes);
}

std::vector<edgeweight>
AffectedNodes::bfsWithoutStartNeighbor(node source, node startNeighbor) {
  edgeweight infDist = std::numeric_limits<edgeweight>::max();

  count z = G.upperNodeIdBound();
  std::vector<edgeweight> distances(G.upperNodeIdBound(), infDist);
  std::vector<bool> visited(z, false);
  std::queue<node> q;

  visited[source] = true;
  distances[source] = 0;

  // Do first iteration separately so we don't have to check every time
  // if we use the edge to be ignored.
  G.forNeighborsOf(source, [&](node v) {
    if (v == startNeighbor || visited[v]) {
      return;
    }
    q.push(v);
    visited[v] = true;
    distances[v] = 1;
  });

  while (!q.empty()) {
    node u = q.front();
    q.pop();
    // insert untouched neighbors into queue
    G.forNeighborsOf(u, [&](node v) {
      if (!visited[v]) {
        q.push(v);
        visited[v] = true;
        distances[v] = distances[u] + 1;
      }
    });
  }

  return distances;
}

std::vector<edgeweight>
AffectedNodes::reverseBfsWithoutStartNeighbor(node source, node startNeighbor) {
  edgeweight infDist = std::numeric_limits<edgeweight>::max();

  count z = G.upperNodeIdBound();
  std::vector<edgeweight> distances(G.upperNodeIdBound(), infDist);
  std::vector<bool> visited(z, false);
  std::queue<node> q;

  visited[source] = true;
  distances[source] = 0;

  // Do first iteration separately so we don't have to check every time
  // if we use the edge to be ignored.
  G.forInNeighborsOf(source, [&](node v) {
    if (v == startNeighbor || visited[v]) {
      return;
    }
    q.push(v);
    visited[v] = true;
    distances[v] = 1;
  });

  while (!q.empty()) {
    node u = q.front();
    q.pop();

    // insert untouched neighbors into queue
    G.forInNeighborsOf(u, [&](node v) {
      if (!visited[v]) {
        q.push(v);
        visited[v] = true;
        distances[v] = distances[u] + 1;
      }
    });
  }

  return distances;
}

std::pair<std::vector<node>, std::vector<edgeweight>>
AffectedNodes::getAffectedNodes(node source,
                                std::vector<edgeweight> &oldDistances,
                                node additionalStartNeighbor) {

  // Run a BFS from the source node but prunes the BFS if the distance to the
  // node has not decreased
  std::queue<node> Q;
  std::vector<edgeweight> newDistances(G.upperNodeIdBound(),
                                       std::numeric_limits<edgeweight>::max());
  std::vector<bool> visited(G.upperNodeIdBound(), false);
  std::vector<node> affectedNodes;

  newDistances[source] = 0;
  visited[source] = true;

  Q.push(source);

  if (additionalStartNeighbor != none) {
    Q.push(additionalStartNeighbor);
    newDistances[additionalStartNeighbor] = 1;
    visited[additionalStartNeighbor] = true;
    affectedNodes.push_back(additionalStartNeighbor);
  }

  do {
    node u = Q.front();
    Q.pop();

    G.forNeighborsOf(u, [&](node v) {
      if (!visited[v]) {
        visited[v] = true;
        newDistances[v] = newDistances[u] + 1;

        if (newDistances[v] < oldDistances[v]) {
          Q.push(v);
          affectedNodes.push_back(v);
        }
      }
    });

  } while (!Q.empty());

  return std::make_pair(affectedNodes, newDistances);
}

std::pair<std::vector<node>, std::vector<edgeweight>>
AffectedNodes::getAffectedNodesBackwards(node source,
                                         std::vector<edgeweight> &oldDistances,
                                         node additionalStartNeighbor) {

  // Run a BFS from the source node but prunes the BFS if the distance to the
  // node has not decreased
  std::queue<node> Q;
  std::vector<edgeweight> newDistances(G.upperNodeIdBound(),
                                       std::numeric_limits<edgeweight>::max());
  std::vector<bool> visited(G.upperNodeIdBound(), false);
  std::vector<node> affectedNodes;

  newDistances[source] = 0;
  visited[source] = true;

  Q.push(source);

  if (additionalStartNeighbor != none) {
    Q.push(additionalStartNeighbor);
    newDistances[additionalStartNeighbor] = 1;
    visited[additionalStartNeighbor] = true;
    affectedNodes.push_back(additionalStartNeighbor);
  }

  do {
    node u = Q.front();
    Q.pop();

    G.forInNeighborsOf(u, [&](node v) {
      if (!visited[v]) {
        visited[v] = true;
        newDistances[v] = newDistances[u] + 1;

        if (newDistances[v] < oldDistances[v]) {
          Q.push(v);
          affectedNodes.push_back(v);
        }
      }
    });

  } while (!Q.empty());

  return std::make_pair(affectedNodes, newDistances);
}

std::vector<edgeweight> AffectedNodes::getDistances() {
  assureFinished();
  return distances;
}

std::vector<node> AffectedNodes::getNodes() {
  assureFinished();
  return nodes;
}

std::vector<edgeweight> AffectedNodes::getImprovements() {
  assureFinished();
  return improvements;
}

} // namespace NetworKit
