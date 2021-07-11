// no-networkit-format
/*
 * SSSP.cpp
 *
 *  Created on: 15.04.2014
 *      Author: cls
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/distance/SSSP.hpp>

namespace NetworKit {

SSSP::SSSP(const Graph &G, node source, bool storePaths,
           bool storeNodesSortedByDistance, node target)
    : Algorithm(), G(&G), source(source), target(target),
      storePaths(storePaths),
      storeNodesSortedByDistance(storeNodesSortedByDistance) {}

const std::vector<edgeweight> &SSSP::getDistances() {
    return distances;
}

std::vector<node> SSSP::getPath(node t, bool forward) const {
    if (!storePaths) {
        throw std::runtime_error("paths have not been stored");
    }
    std::vector<node> path;
    if (previous[t].empty()) { // t is not reachable from source
        WARN("there is no path from ", source, " to ", t);
        return path;
    }
    node v = t;
    while (v != source) {
        path.push_back(v);
        v = previous[v].front();
    }
    path.push_back(source);

    if (forward) {
        std::reverse(path.begin(), path.end());
    }
    return path;
}

std::set<std::vector<node>> SSSP::getPaths(node t, bool forward) const {

    std::set<std::vector<node>> paths;
    if (previous[t].empty()) { // t is not reachable from source
        WARN("there is no path from ", source, " to ", t);
        return paths;
    }

    std::vector<node> targetPredecessors = previous[t];

#pragma omp parallel for schedule(dynamic)
    for (omp_index i = 0; i < static_cast<omp_index>(targetPredecessors.size());
         ++i) {

        std::stack<std::vector<node>> stack;
        std::vector<std::vector<node>> currPaths;

        node pred = targetPredecessors[i];
        if (pred == source) {
            currPaths.push_back({t, pred});
        } else {
            stack.push({t, pred});
        }

        while (!stack.empty()) {

            node topPathLastNode = stack.top().back();

            if (topPathLastNode == source) {
                currPaths.push_back(stack.top());
                stack.pop();
                continue;
            }

            std::vector<node> topPath = stack.top();
            stack.pop();

            std::vector<node> currPredecessors = previous[topPath.back()];

            for (node currPredecessor : currPredecessors) {
                std::vector<node> suffix(topPath);
                suffix.push_back(currPredecessor);
                stack.push(suffix);
            }
        }

#pragma omp critical
        paths.insert(currPaths.begin(), currPaths.end());
    }

    if (forward) {
        std::set<std::vector<node>> reversedPaths;
        for (std::vector<node> path : paths) {
            std::reverse(std::begin(path), std::end(path));
            reversedPaths.insert(path);
        }
        paths = reversedPaths;
    }

    return paths;
}

const std::vector<node> &SSSP::getNodesSortedByDistance() const {
    if (!storeNodesSortedByDistance) {
        throw std::runtime_error(
            "Nodes sorted by distance have not been stored. Set "
            "storeNodesSortedByDistance in the constructor to true to enable "
            "this behaviour.");
    }

    assert(!nodesSortedByDistance.empty());
    return nodesSortedByDistance;
}

} /* namespace NetworKit */
