/*
 * ApproximatePageRank.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include <queue>

#include <networkit/scd/ApproximatePageRank.hpp>

namespace NetworKit {

ApproximatePageRank::ApproximatePageRank(const Graph &g, double alpha, double epsilon)
    : g(&g), alpha(alpha), eps(epsilon) {}

std::vector<std::pair<node, double>> ApproximatePageRank::run(const std::set<node> &seeds) {
    double initRes = 1.0 / seeds.size();
    std::queue<node> activeNodes;
    for (node s : seeds) {
        prRes[s] = std::make_pair(0.0, initRes);
        activeNodes.push(s);
    }

    auto push = [&](const node u, std::queue<node> &activeNodes) {
        double res = prRes[u].second;
        double volume = g->weightedDegree(u, true);

        g->forNeighborsOf(u, [&](node, const node v, const edgeweight w) {
            double mass = (1.0 - alpha) * res * w / (2.0 * volume);
            double volV = g->weightedDegree(v, true);
            // the first check is for making sure the node is not added twice.
            // the second check ensures that enough residual is left.
            if (prRes[v].second < volV * eps && (prRes[v].second + mass) >= eps * volV) {
                activeNodes.push(v);
            }
            prRes[v].second += mass;
        });

        prRes[u] = std::make_pair(prRes[u].first + alpha * res, (1.0 - alpha) * res / 2);
        if ((prRes[u].second / volume) >= eps) {
            activeNodes.push(u);
        }
    };

    while (!activeNodes.empty()) {
        node v = activeNodes.front();
        activeNodes.pop();
        push(v, activeNodes);
    }

    std::vector<std::pair<node, double>> pr;
    pr.reserve(prRes.size());

    for (auto it = prRes.begin(); it != prRes.end(); it++) {
        pr.emplace_back(it->first, it->second.first);
    }

    return pr;
}

std::vector<std::pair<node, double>> ApproximatePageRank::run(node seed) {
    return run(std::set<node>{seed});
}

} /* namespace NetworKit */
