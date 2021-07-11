// no-networkit-format
/*
 * BiconnectedeComponents.cpp
 *
 * Created on: March 2018
 * 		 Author: Eugenio Angriman
 */

#include <cmath>
#include <stack>

#include <networkit/components/BiconnectedComponents.hpp>

namespace NetworKit {

BiconnectedComponents::BiconnectedComponents(const Graph &G) : G(&G) {
  if (G.isDirected()) {
    throw std::runtime_error(
        "Error, biconnected components cannot be computed on directed graphs.");
  }
}

void BiconnectedComponents::init() {
  n = G->numberOfNodes();
  idx = 0;
  nComp = 0;
  level.assign(n, 0);
  lowpt.assign(n, 0);
  parent.assign(n, 0);
  visited.assign(n, false);
  isRoot.assign(n, false);
  componentsOfNode.clear();
  componentsOfNode.resize(n);
}

void BiconnectedComponents::run() {

  init();

  auto visitNode = [&](node u) {
    level[u] = idx;
    lowpt[u] = idx;
    ++idx;
    visited[u] = true;
  };

  std::stack<std::pair<node, Graph::NeighborIterator>> stack;
  std::vector<std::pair<node, node>> edgeStack;
  G->forNodes([&](node v) {
    if (visited[v]) {
      return;
    }

    isRoot[v] = true;
    stack.push(std::make_pair(v, G->neighborRange(v).begin()));

    do {
      node u = stack.top().first;
      auto &iter = stack.top().second;
      if (!visited[u]) {
        visitNode(u);
      }

      for (; iter != G->neighborRange(u).end(); ++iter) {
        node neighbor = *iter;
        if (!visited[neighbor]) {
          visitNode(neighbor);
          parent[neighbor] = u;
          stack.push(std::make_pair(neighbor, G->neighborRange(neighbor).begin()));
          edgeStack.emplace_back(u, neighbor);
          break;
        } else if (neighbor != parent[u] && level[neighbor] < level[u]) {
            edgeStack.emplace_back(u, neighbor);
            lowpt[u] = std::min(lowpt[u], level[neighbor]);
        }
      }

      if (iter == G->neighborRange(u).end()) {
        stack.pop();

        if (isRoot[u]) {
          continue;
        }

        node v = parent[u];
        lowpt[v] = std::min(lowpt[v], lowpt[u]);
        if (lowpt[u] >= level[v]) {
          auto top = edgeStack.back();
          componentSizes.insert(std::make_pair(nComp, 0));
          while (level[top.first] >= level[u]) {
            edgeStack.pop_back();
            newComponent(top);
            top = edgeStack.back();
          }

          for (auto rit = edgeStack.rbegin(); rit != edgeStack.rend(); ++rit) {
            if ((*rit).first == v && (*rit).second == u) {
              newComponent(*rit);
              edgeStack.erase(std::next(rit).base());
              break;
            }
          }
          ++nComp;
        }
      }
    } while (!stack.empty());

    edgeStack.clear();
  });

  hasRun = true;
}

void BiconnectedComponents::newComponent(std::pair<node, node> e) {
  node u = e.first;
  node v = e.second;
  count prevSize = componentsOfNode[u].size() + componentsOfNode[v].size();
  componentsOfNode[u].insert(nComp);
  componentsOfNode[v].insert(nComp);
  componentSizes[nComp] +=
      componentsOfNode[u].size() + componentsOfNode[v].size() - prevSize;
}
} // namespace NetworKit
