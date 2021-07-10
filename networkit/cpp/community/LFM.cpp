/*
 * LFM.cpp
 *
 *  Created on: 10.12.2020
 *      Author: John Gelhausen
 */
// networkit-format

#include <iostream>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/community/LFM.hpp>

namespace NetworKit {

LFM::LFM(const Graph &G, SelectiveCommunityDetector &scd)
    : OverlappingCommunityDetectionAlgorithm(G), scd(&scd) {}

void LFM::run() {
    Aux::SignalHandler handler;

    Cover zeta(G->upperNodeIdBound());
    index o = 0;
    zeta.setUpperBound(o);

    G->forNodesInRandomOrder([&](node u) {
        handler.assureRunning();
        if (!zeta.contains(u)) {
            std::set<node> community = scd->expandOneCommunity(u);
            o++;
            zeta.setUpperBound(o);

            handler.assureRunning();

            for (node n : community) {
                zeta.addToSubset(o - 1, n);
            }
        }
    });

    result = std::move(zeta);
    hasRun = true;
}

} /* namespace NetworKit */
