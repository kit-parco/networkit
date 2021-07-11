// no-networkit-format
/*
 * APSP.cpp
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe
 */

#include <omp.h>

#include <networkit/distance/APSP.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Dijkstra.hpp>

namespace NetworKit {

APSP::APSP(const Graph &G) : Algorithm(), G(G) {}

void APSP::run() {
    count n = G.upperNodeIdBound();
    std::vector<edgeweight> distanceVector(n, 0.0);
    distances.resize(n, distanceVector);

    count nThreads = omp_get_max_threads();
    sssps.resize(nThreads);
#pragma omp parallel
    {
        omp_index i = omp_get_thread_num();
        if (G.isWeighted())
            sssps[i] = std::unique_ptr<SSSP>(new Dijkstra(G, 0, false));
        else
            sssps[i] = std::unique_ptr<SSSP>(new BFS(G, 0, false));
    }

#pragma omp parallel for schedule(dynamic)
    for (omp_index source = 0; source < n; ++source) {
        auto sssp = sssps[omp_get_thread_num()].get();
        sssp->setSource(source);
        sssp->run();
        distances[source] = sssp->getDistances();
    }

    hasRun = true;
}

std::string APSP::toString() const {
    return "All-Pairs Shortest Path Algorithm";
}

} /* namespace NetworKit */
