// no-networkit-format
/*
 * LaplacianCentrality.cpp
 *
 *  Created on: 08.03.2018
 *      Author: Kolja Esders
 */

#include <networkit/centrality/LaplacianCentrality.hpp>

namespace NetworKit {

LaplacianCentrality::LaplacianCentrality(const Graph& G, bool normalized): Centrality(G, normalized) {
}

void LaplacianCentrality::run() {
    scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);
    double totalLaplacianEnergy = 0.0;

    G.parallelForNodes([&](node u) {
        count degreeU = G.weightedDegree(u);
        double energyLossOnNodeDrop = static_cast<double>(degreeU * degreeU);
#pragma omp atomic
        totalLaplacianEnergy += static_cast<double>(degreeU * degreeU);

        G.forNeighborsOf(u, [&](node v, edgeweight ew) {
            energyLossOnNodeDrop += ew * (ew + 2 * G.weightedDegree(v));
#pragma omp atomic
            totalLaplacianEnergy += ew * ew;
        });

        scoreData[u] = energyLossOnNodeDrop;
    });

    if (!normalized) {
        hasRun = true;
        return;
    }

    G.parallelForNodes([&](node u) {
        scoreData[u] /= totalLaplacianEnergy;
    });

    hasRun = true;
}

} /* namespace NetworKit */
