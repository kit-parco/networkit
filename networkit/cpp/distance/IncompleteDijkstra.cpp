/*
 * IncompleteDijkstra.cpp
 *
 *  Created on: 15.07.2014
 *      Author: dhoske
 */

#include <limits>
#include <type_traits>

#include <networkit/distance/IncompleteDijkstra.hpp>

namespace NetworKit {

IncompleteDijkstra::IncompleteDijkstra(const Graph *G, const std::vector<node> &sources,
                                       const std::unordered_set<node> *explored)
    : G(G), explored(explored), heap(CompareDistance(&dists)) {
    if (!G) {
        throw std::invalid_argument("G is null");
    }

    dists.resize(G->upperNodeIdBound(), std::numeric_limits<edgeweight>::max());
    heap.reserve(G->upperNodeIdBound());

    for (const auto source : sources) {
        if (!explored || explored->find(source) == explored->end()) {
            dists[source] = 0.0;
            heap.update(source);
        }
    }
}

bool IncompleteDijkstra::hasNext() {
    return !heap.empty();
}

std::pair<node, edgeweight> IncompleteDijkstra::next() {
    if (!hasNext()) {
        throw std::invalid_argument("No next element");
    }

    // Extract nearest node
    const auto u = heap.extract_top();
    const auto dist_u = dists[u];

    // Relax all of its edges
    G->forNeighborsOf(u, [&](const node v, const edgeweight dist_uv) {
        if (explored && explored->find(v) != explored->end()) {
            return;
        }

        if (dist_u + dist_uv < dists[v]) {
            dists[v] = dist_u + dist_uv;
            heap.update(v);
        }
    });

    return {u, dist_u};
}

} // namespace NetworKit
