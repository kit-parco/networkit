/*
 * AllSimplePaths.cpp
 *
 *  Created on: 23.06.2017
 *      Author: Eugenio Angriman
 */

#include <memory>
#include <omp.h>

#include <networkit/reachability/AllSimplePaths.hpp>

namespace NetworKit {

AllSimplePaths::AllSimplePaths(const Graph &G, node source, node target, count cutoff)
    : G(&G), source(source), target(target), cutoff(cutoff) {
    if (!G.isDirected())
        throw std::runtime_error(
            "Error, AllSimplePaths class has been implemented for directed graphs only.");

    if (G.isWeighted())
        throw std::runtime_error(
            "Error, AllSimplePaths class has been implemented for unweighted graphs only.");

    if (!G.hasNode(source))
        throw std::runtime_error("Error, source node not in graph.");

    if (!G.hasNode(target))
        throw std::runtime_error("Error, source node not in graph.");

    if (source == target)
        throw std::runtime_error("Error, source is equal to the target.");

    if (cutoff < 1)
        throw std::runtime_error("Error, cutoff = 0.");

    distanceToTarget.assign(G.upperNodeIdBound(), none);
    distanceFromSource.assign(G.upperNodeIdBound(), none);
}

void AllSimplePaths::run() {

    std::queue<node> q;
    q.push(target);
    distanceToTarget[target] = 0;
    distanceFromSource[source] = 0;
    bool reachable = false;

    // Reverse BFS from target. Labels all nodes with their distance to target
    do {
        const node curr = q.front();
        q.pop();

        // Distance neighbor nodes to target
        const count d = distanceToTarget[curr] + 1;

        // Labels all neighbors with their distances to the target
        // only if their distance to target is inside cutoff.
        if (d <= cutoff) {

            // Reverse visit with InNeighbors.
            G->forInNeighborsOf(curr, [&](node v) {
                // Source node reaches target node.
                if (v == source) {
                    if (d > cutoff)
                        throw std::runtime_error(
                            "Error, source node cannot reach target node within the given cutoff.");

                    reachable = true;
                }
                // Unexplored node.
                if (distanceToTarget[v] == none) {
                    distanceToTarget[v] = d;
                    q.push(v);
                }
            });
        }

    } while (!q.empty());

    // Source cannot reach target. Stopping algorithm.
    if (!reachable)
        throw std::runtime_error("Error, source node cannot reach target node");

    // Start another BFS from source.
    q.push(source);

    // BFS from source. Labels all nodes with their distance from source.
    do {
        const node curr = q.front();
        q.pop();

        // Distance neighbor nodes from source
        count d = distanceFromSource[curr] + 1;

        // Labels all neighbors with their distances from source
        // only if they can be part of a path from source to target.
        G->forNeighborsOf(curr, [&](node v) {
            if (distanceFromSource[v] == none && distanceToTarget[v] != none
                && (cutoff == none || d + distanceToTarget[v] <= cutoff)) {
                distanceFromSource[v] = d;
                q.push(v);
            }
        });
    } while (!q.empty());

    // Computing all simple paths.
    computePaths();

    hasRun = true;
}

void AllSimplePaths::computePaths() {

    std::vector<std::vector<node>> availableSources(G->upperNodeIdBound());
    G->parallelForNodes([&](node v) {
        if (v != target && (cutoff == none || distanceFromSource[v] != none))
            availableSources[v] = getAvailableSources(v);
    });

    omp_lock_t lock;
    omp_init_lock(&lock);
    paths.reserve(availableSources[source].size());

    using StackElem = std::pair<std::vector<node>, std::vector<bool>>;
#pragma omp parallel for schedule(dynamic)
    for (omp_index i = 0; i < static_cast<omp_index>(availableSources[source].size()); ++i) {
        std::vector<StackElem> stack;
        std::vector<bool> visited(G->upperNodeIdBound(), false);
        visited[source] = true;
        {
            std::vector<node> v{{source, availableSources[source][i]}};
            stack.emplace_back(v, visited);
        }

        std::vector<std::vector<node>> currPaths;

        do {
            const node curr = stack.back().first.back();

            if (curr == target) {
                currPaths.push_back(stack.back().first);
                stack.pop_back();
                continue;
            }

            // Pointer to last vector inserted in the stack.
            StackElem *currPair = &stack.back();

            // Erase paths that reached the cutoff without reaching target.
            if (cutoff != none && currPair->first.size() - 1 >= cutoff) {
                stack.pop_back();
                continue;
            }

            count prevSize = stack.size();
            currPair->second[curr] = true;
            bool toEnqueue = false;

            if (!availableSources[curr].empty()) {
                node s = availableSources[curr][0];
                if (!currPair->second[s])
                    currPair->first.push_back(s);
                else
                    toEnqueue = true;
            }

            count step = 0;
            stack.reserve(availableSources.size() - 1);
            currPair = &stack.back();

            for (count i = 1; i < availableSources[curr].size(); ++i) {

                const node s = availableSources[curr][i];
                if (!currPair->second[s]) { // Not visited. Avoid loops.
                    // Not in stack. Creating new pair.
                    if (toEnqueue) {
                        currPair->first.push_back(s);
                        toEnqueue = false;
                    } else {
                        stack.emplace_back(*currPair);
                        auto &newPair = stack.back();
                        newPair.first[newPair.first.size() - 1] = s;
                        ++step;
                        currPair = &stack[stack.size() - 1 - step];
                    }
                }
            }

            // No new available sources found for this prefix.
            if (toEnqueue && prevSize == stack.size())
                stack.pop_back();

        } while (!stack.empty());

        omp_set_lock(&lock);
        paths.insert(paths.end(), currPaths.begin(), currPaths.end());
        omp_unset_lock(&lock);
    }
}

std::vector<node> AllSimplePaths::getAvailableSources(node s, count pathLength) {
    std::vector<node> availableSources;
    G->forNeighborsOf(s, [&](node v) {
        // Make sure we are visiting a node that can reach target. Avoid to consider source as a
        // possible source.
        if (distanceFromSource[v] != none)
            if (cutoff == none || pathLength + distanceToTarget[v] <= cutoff)
                availableSources.push_back(v);
    });

    return availableSources;
}

} /* namespace NetworKit */
