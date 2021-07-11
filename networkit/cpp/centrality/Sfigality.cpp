// no-networkit-format
/*
 * Sfigality.cpp
 *
 *  Created on: 20.01.2016
 *      Author: Christian Staudt
 */

#include <networkit/centrality/Sfigality.hpp>
#include <algorithm>

namespace NetworKit {

Sfigality::Sfigality(const Graph& G) : Centrality(G, true) {
}

void Sfigality::run() {
    scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);

    G.parallelForNodes([&](node u) {
        count sf = 0;
        G.forEdgesOf(u, [&](node u, node v){
            if (G.degree(u) < G.degree(v)) {
                sf += 1;
            };
        });
        scoreData[u] = sf / (double) G.degree(u);
    });

    hasRun = true;
}


double Sfigality::maximum() {
    throw std::runtime_error("Not Implemented");
    return G.isEmpty() ? 0. : static_cast<double>(G.numberOfNodes() - 1);
}


} /* namespace NetworKit */
