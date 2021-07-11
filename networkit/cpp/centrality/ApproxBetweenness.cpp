// no-networkit-format
/*
 * ApproxBetweenness.cpp
 *
 *  Created on: 09.04.2014
 *      Author: cls
 */

#include <networkit/centrality/ApproxBetweenness.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/SSSP.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <omp.h>

namespace NetworKit {

ApproxBetweenness::ApproxBetweenness(const Graph& G, double epsilon, double delta, double universalConstant) : Centrality(G, true), epsilon(epsilon), delta(delta), universalConstant(universalConstant) {

}


void ApproxBetweenness::run() {
    Aux::SignalHandler handler;
    scoreData.clear();
    scoreData.resize(G.upperNodeIdBound());

    edgeweight vd = 0;

    Diameter diam(G, DiameterAlgo::estimatedPedantic);
    diam.run();
    vd = static_cast<edgeweight>(diam.getDiameter().first);

    if (vd <= 2) {
        hasRun = true;
        return;
    }

    INFO("estimated diameter: ", vd);
    r = ceil((universalConstant / (epsilon * epsilon)) * (floor(log2(vd - 2)) + 1 - log(delta)));

    INFO("taking ", r, " path samples");
    handler.assureRunning();
    #pragma omp parallel
    {
        auto ssspPtr = G.isWeighted() ? std::unique_ptr<SSSP>(new Dijkstra(G, 0, true, false))
                                      : std::unique_ptr<SSSP>(new BFS(G, 0, true, false));

#pragma omp for
        for (omp_index i = 1; i <= static_cast<omp_index>(r); i++) {
            DEBUG("sample ", i);
            // sample random node pair
            node u, v;
            u = GraphTools::randomNode(G);
            do {
                v = GraphTools::randomNode(G);
            } while (v == u);

            auto &sssp = *ssspPtr;
            sssp.setSource(u);
            sssp.setTarget(v);

            DEBUG("running shortest path algorithm for node ", u);
            if (!handler.isRunning()) continue;
            sssp.run();
            if (!handler.isRunning()) continue;
            if (sssp.numberOfPaths(v) > 0) { // at least one path between {u, v} exists
                DEBUG("updating estimate for path ", u, " <-> ", v);
                // random path sampling and estimation update
                node t = v;
                while (t != u)  {
                    // sample z in P_u(t) with probability sigma_uz / sigma_us
                    std::vector<std::pair<node, double> > choices;
                    for (node z : sssp.getPredecessors(t)) {
                        bigfloat tmp = sssp.numberOfPaths(z) / sssp.numberOfPaths(t);
                        double weight;
                        tmp.ToDouble(weight);
                        choices.emplace_back(z, weight); 	// sigma_uz / sigma_us
                    }
                    node z = Aux::Random::weightedChoice(choices);
                    assert (z <= G.upperNodeIdBound());
                    if (z != u)
#pragma omp atomic
                        scoreData[z] += 1. / static_cast<double>(r);
                    t = z;
                }
            }
        }
    }
    handler.assureRunning();

    hasRun = true;
}


count ApproxBetweenness::numberOfSamples() const {
    assureFinished();
    INFO("Estimated number of samples", r);
    return r;
}


} /* namespace NetworKit */
