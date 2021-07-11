
#ifndef NETWORKIT_GRAPH_BFS_HPP_
#define NETWORKIT_GRAPH_BFS_HPP_

#include <array>
#include <queue>
#include <vector>

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

namespace Traversal {

/**
 * Calls the given BFS handle with distance parameter
 */
template <class F>
auto callBFSHandle(F &f, node u, count dist) -> decltype(f(u, dist)) {
    return f(u, dist);
}

/**
 * Calls the given BFS handle without distance parameter
 */
template <class F>
auto callBFSHandle(F &f, node u, count) -> decltype(f(u)) {
    return f(u);
}

/**
 * Iterate over nodes in breadth-first search order starting from the nodes within the given range.
 *
 * @param G The input graph.
 * @param first The first element of the range.
 * @param last The end of the range.
 * @param handle Takes a node as input parameter.
 */
template <class InputIt, typename L>
void BFSfrom(const Graph &G, InputIt first, InputIt last, L handle) {
    std::vector<bool> marked(G.upperNodeIdBound());
    std::queue<node> q, qNext;
    count dist = 0;
    // enqueue start nodes
    for (; first != last; ++first) {
        q.push(*first);
        marked[*first] = true;
    }
    do {
        const auto u = q.front();
        q.pop();
        // apply function
        callBFSHandle(handle, u, dist);
        G.forNeighborsOf(u, [&](node v) {
            if (!marked[v]) {
                qNext.push(v);
                marked[v] = true;
            }
        });
        if (q.empty() && !qNext.empty()) {
            q.swap(qNext);
            ++dist;
        }
    } while (!q.empty());
}

/**
 * Iterate over nodes in breadth-first search order starting from the given source node.
 *
 * @param G The input graph.
 * @param source The source node.
 * @param handle Takes a node as input parameter.
 */
template <typename L>
void BFSfrom(const Graph &G, node source, L handle) {
    std::array<node, 1> startNodes{{source}};
    BFSfrom(G, startNodes.begin(), startNodes.end(), handle);
}

/**
 * Iterate over edges in breadth-first search order starting from the given source node.
 *
 * @param G The input graph.
 * @param source The source node.
 * @param handle Takes a node as input parameter.
 */
template <typename L>
void BFSEdgesFrom(const Graph &G, node source, L handle) {
    std::vector<bool> marked(G.upperNodeIdBound());
    std::queue<node> q;
    q.push(source); // enqueue root
    marked[source] = true;
    do {
        const auto u = q.front();
        q.pop();
        // apply function
        G.forNeighborsOf(u, [&](node, node v, edgeweight w, edgeid eid) {
            if (!marked[v]) {
                handle(u, v, w, eid);
                q.push(v);
                marked[v] = true;
            }
        });
    } while (!q.empty());
}

} // namespace Traversal

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_BFS_HPP_
