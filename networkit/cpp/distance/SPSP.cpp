/*
 * SPSP.cpp
 *
 *  Created on: 23.10.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#include <memory>
#include <omp.h>

#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/SPSP.hpp>
#include <networkit/distance/SSSP.hpp>

namespace NetworKit {

void SPSP::run() {
    distances.resize(sources.size());

#pragma omp parallel
    {
        std::unique_ptr<SSSP> sssp;
        if (G->isWeighted())
            sssp = std::unique_ptr<SSSP>(new Dijkstra(*G, 0, false));
        else
            sssp = std::unique_ptr<SSSP>(new BFS(*G, 0, false));

#pragma omp for
        for (omp_index i = 0; i < sources.size(); ++i) {
            sssp->setSource(sources[i]);
            sssp->run();
            distances[i] = sssp->getDistances();
        }
    }

    hasRun = true;
}

} // namespace NetworKit
