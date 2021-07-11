/*
 * BidirectionalDijkstra.cpp
 *
 *  Created on: 14.06.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#include <limits>

#include <networkit/distance/BidirectionalDijkstra.hpp>

namespace NetworKit {

void BidirectionalDijkstra::run() {
    if (!G->isWeighted())
        WARN("The graph is unweighted, use BidirectionalBFS for better "
             "efficiency!");
    stDist = 0.;
    if (source == target) {
        hasRun = true;
        return;
    }

    init();
    visited.resize(G->upperNodeIdBound(), ts);
    if (ts++ == 128) {
        ts = 1;
        std::fill(visited.begin(), visited.end(), 0);
    }
    count n = G->upperNodeIdBound();
    if (storePred) {
        pred[source] = source;
        predT.resize(n);
        predT[target] = target;
    }

    edgeweight infdist = std::numeric_limits<edgeweight>::max();

    std::fill(dist1.begin(), dist1.end(), infdist);
    dist1.resize(n, infdist);
    std::fill(dist2.begin(), dist2.end(), infdist);
    dist2.resize(n, infdist);

    dist1[source] = dist2[target] = 0.;
    // Necessary only for the target if the target ball is never expanded.
    visited[target] = ts + ballMask;

    h1.clear();
    h1.reserve(n);
    h1.push(source);

    h2.clear();
    h2.reserve(n);
    h2.push(target);

    stH.clear();
    stH.reserve(n);
    std::stack<std::deque<node>> pathsToBuild;

    // Expands a ball; idx is the ball index (either 0 or 1<<7)
    auto expand = [&](tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<edgeweight>> &h,
                      uint8_t idx) {
        node u = h.extract_top();
        if (ts != (visited[u] & ~ballMask)) {
            // u not visited yet
            visited[u] = ts + idx;
        } else {
            // u is in the other ball
            if ((visited[u] & ballMask) != idx)
                return u;
        }

        auto visitEdge = [&](node v, edgeweight w) {
            edgeweight du = idx ? dist2[u] : dist1[u];
            edgeweight &dv = idx ? dist2[v] : dist1[v];
            if (dv > du + w) {
                dv = du + w;
                h.update(v);
                (idx ? predT : pred)[v] = u;
                // v is in the other heap
                if ((idx ? h1 : h2).contains(v))
                    stH.update(v);
            }
        };

        if (!G->isDirected() || idx != ballMask)
            G->forNeighborsOf(u, visitEdge);
        else
            G->forInNeighborsOf(u, visitEdge);

        return none;
    };

    node u = none;

    do {
        if (h1.size() <= h2.size())
            u = expand(h1, 0);
        else
            u = expand(h2, ballMask);
    } while (u == none && !h1.empty() && !h2.empty());

    if (u != none) {
        stDist = dist1[u] + dist2[u];
        auto sumDists = [&](node x) { return dist1[x] + dist2[x]; };

        if (!stH.empty() && stDist > sumDists(stH.top())) {
            stDist = sumDists(stH.top());
            u = stH.extract_top();
        }

        if (storePred) {
            while (u != target) {
                node u1 = predT[u];
                pred[u1] = u;
                u = u1;
            }
            buildPath();
        }
    } else {
        stDist = infdist;
        WARN("Source cannot reach target!");
    }

    hasRun = true;
}
} // namespace NetworKit
