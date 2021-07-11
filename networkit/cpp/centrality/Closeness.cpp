// no-networkit-format
/*
 * Closeness.cpp
 *
 *  Created on: 03.10.2014
 *     Authors: nemes
 *              Eugenio Angriman <angrimae@hu-berlin.de>
 */

#include <omp.h>
#include <queue>

#include <networkit/centrality/Closeness.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>

namespace NetworKit {

Closeness::Closeness(const Graph &G, bool normalized,
                     const ClosenessVariant variant)
    : Centrality(G, normalized), variant(variant) {
    if (variant == ClosenessVariant::standard) {
        checkConnectedComponents();
    }
}

Closeness::Closeness(const Graph &G, bool normalized, bool checkConnectedness)
    : Centrality(G, normalized), variant(ClosenessVariant::standard) {
    if (checkConnectedness)
        checkConnectedComponents();
}

void Closeness::checkConnectedComponents() const {
    bool multipleComponents = false;
    if (G.isDirected()) {
        StronglyConnectedComponents scc(G);
        scc.run();
        multipleComponents = scc.numberOfComponents() > 1;
    } else {
        ConnectedComponents cc(G);
        cc.run();
        multipleComponents = cc.numberOfComponents() > 1;
    }
    if (multipleComponents) {
        throw std::runtime_error(
            "Error: the standard definition of closeness is not defined on "
            "disconnected graphs. On disconnected graphs, use the generalized "
            "definition instead.");
    }
}

void Closeness::run() {
    count n = G.upperNodeIdBound();

    scoreData.clear();
    scoreData.resize(n);
    visited.clear();
    visited.resize(omp_get_max_threads(), std::vector<uint8_t>(n));
    ts.clear();
    ts.resize(omp_get_max_threads(), 0);

    if (G.isWeighted()) {
        dDist.resize(omp_get_max_threads(), std::vector<double>(n));
        heaps.reserve(omp_get_max_threads());
        for (int i = 0; i < omp_get_max_threads(); ++i) {
            heaps.emplace_back((Compare(dDist[i])));
            heaps.back().reserve(n);
        }
        dijkstra();
    } else {
        uDist.resize(omp_get_max_threads(), std::vector<count>(n));
        bfs();
    }

    hasRun = true;
}

void Closeness::bfs() {
#pragma omp parallel for schedule(dynamic)
    for (omp_index u = 0; u < static_cast<omp_index>(G.upperNodeIdBound());
         ++u) {
        uint8_t &ts_ = ts[omp_get_thread_num()];
        auto &dist_ = uDist[omp_get_thread_num()];
        auto &visited_ = visited[omp_get_thread_num()];

        if (ts_++ == 255 && variant != ClosenessVariant::standard) {
            ts_ = 1;
            std::fill(visited_.begin(), visited_.end(), 0);
        }

        std::queue<node> q;
        q.push(u);
        dist_[u] = 0;
        visited_[u] = ts_;
        double sum = 0.;
        count reached = 1;

        do {
            node x = q.front();
            q.pop();
            G.forNeighborsOf(x, [&](node y) {
                if (visited_[y] != ts_) {
                    visited_[y] = ts_;
                    dist_[y] = dist_[x] + 1;
                    sum += dist_[y];
                    ++reached;
                    q.push(y);
                }
            });

        } while (!q.empty());

        updateScoreData(u, reached, sum);
    }
}

void Closeness::dijkstra() {
#pragma omp parallel for schedule(dynamic)
    for (omp_index u = 0; u < static_cast<omp_index>(G.upperNodeIdBound());
         ++u) {
        uint8_t &ts_ = ts[omp_get_thread_num()];
        auto &dist_ = dDist[omp_get_thread_num()];
        auto &visited_ = visited[omp_get_thread_num()];
        auto &heap = heaps[omp_get_thread_num()];

        if (ts_++ == 255 && variant != ClosenessVariant::standard) {
            ts_ = 1;
            std::fill(visited_.begin(), visited_.end(), 0);
        }

        heap.push(u);
        dist_[u] = 0.;
        visited_[u] = ts_;

        double sumDist = 0.;
        count reached = 1;

        do {
            const auto x = heap.extract_top();
            sumDist += dist_[x];
            G.forNeighborsOf(x, [&](node y, edgeweight ew) {
                if (ts_ != visited_[y]) {
                    ++reached;
                    visited_[y] = ts_;
                    dist_[y] = dist_[x] + ew;
                    heap.push(y);
                } else if (dist_[y] > dist_[x] + ew) {
                    dist_[y] = dist_[x] + ew;
                    heap.update(y);
                }
            });
        } while (!heap.empty());

        updateScoreData(u, reached, sumDist);
    }
}

} // namespace NetworKit
