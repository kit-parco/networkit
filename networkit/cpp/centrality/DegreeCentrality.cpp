// no-networkit-format
/*
 * DegreeCentrality.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <networkit/centrality/DegreeCentrality.hpp>

namespace NetworKit {

DegreeCentrality::DegreeCentrality(const Graph& G, bool normalized, bool outDeg, bool ignoreSelfLoops) : Centrality(G, normalized), outDeg(outDeg), ignoreSelfLoops(ignoreSelfLoops) {
}

void DegreeCentrality::run() {
    scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);
    ignoreSelfLoops = ignoreSelfLoops && (G.numberOfSelfLoops() > 0); //this way we don't need to check for self loops if there are none
    
    if (G.isDirected() && !outDeg) {
        G.parallelForNodes([&](node u) {
            scoreData[u] = G.degreeIn(u);
            if(ignoreSelfLoops && G.hasEdge(u, u)) scoreData[u] -= 1;
        });
    } else {
        G.parallelForNodes([&](node u) {
            scoreData[u] = G.degree(u);
            if(ignoreSelfLoops && G.hasEdge(u, u)) scoreData[u] -= 1;
        });
    }

    if (normalized) {
        const double maxDeg = maximum();
        G.parallelForNodes([&](node u) {
            scoreData[u] = scoreData[u] / maxDeg;
        });
    }
    hasRun = true;
}


double DegreeCentrality::maximum() {
    if (G.isEmpty())
        return 0.;
    if (ignoreSelfLoops) {
        return static_cast<double>(G.numberOfNodes() - 1);
    } else {
        return static_cast<double>(G.numberOfNodes());
    }
}


} /* namespace NetworKit */
