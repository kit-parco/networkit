
#ifndef NETWORKIT_GRAPH_DIJKSTRA_HPP_
#define NETWORKIT_GRAPH_DIJKSTRA_HPP_

#include <limits>
#include <vector>

#include <networkit/graph/Graph.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

namespace NetworKit {

namespace Traversal {

/**
 * Calls the given Dijkstra handle with distance parameter
 */
template <class F>
auto callDijkstraHandle(F &f, node u, edgeweight dist) -> decltype(f(u, dist)) {
    return f(u, dist);
}

/**
 * Calls the given Dijkstra handle without distance parameter
 */
template <class F>
auto callDijkstraHandle(F &f, node u, edgeweight) -> decltype(f(u)) {
    return f(u);
}

/**
 * Iterate over nodes with Dijkstra starting from the nodes within the given range.
 *
 * @param G The input graph.
 * @param first The first element of the range.
 * @param last The end of the range.
 * @param lambda Takes a node and its distance from the nodes in the range as input parameters.
 */
template <class InputIt, typename Handle>
void DijkstraFrom(const Graph &G, InputIt first, InputIt last, Handle handle) {
    std::vector<edgeweight> distance(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max());
    const auto compareDistance = [&distance](node u, node v) noexcept -> bool {
        return distance[u] < distance[v];
    };

    tlx::d_ary_addressable_int_heap<node, 2, decltype(compareDistance)> heap{compareDistance};

    std::for_each(first, last, [&heap, &distance](const node u) {
        heap.push(u);
        distance[u] = 0;
    });

    while (!heap.empty()) {
        const auto u = heap.extract_top();
        callDijkstraHandle(handle, u, distance[u]);
        G.forNeighborsOf(u, [&](const node v, const edgeweight w) {
            if (distance[v] > distance[u] + w) {
                distance[v] = distance[u] + w;
                heap.update(v);
            }
        });
    }
}

/**
 * Iterate over nodes with Dijkstra starting from the nodes within the given range.
 *
 * @param G The input graph.
 * @param u The source node.
 * @param lambda Takes a node and its distance from `u` as input parameters.
 */
template <typename Lambda>
void DijkstraFrom(const Graph &G, node u, Lambda lambda) {
    std::vector<node> vec({u});
    DijkstraFrom(G, vec.begin(), vec.end(), lambda);
}

} // namespace Traversal

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_DIJKSTRA_HPP_
