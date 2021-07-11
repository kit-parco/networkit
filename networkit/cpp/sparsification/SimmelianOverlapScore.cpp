// no-networkit-format
/*
 * SimmelianOverlapScore.cpp
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#include <networkit/sparsification/SimmelianOverlapScore.hpp>
#include <limits>

namespace NetworKit {

SimmelianOverlapScore::SimmelianOverlapScore(const Graph& G, const std::vector<count>& triangles, count maxRank) :
        SimmelianScore(G, triangles), maxRank(maxRank) {}

void SimmelianOverlapScore::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(*G, *triangles);
    scoreData.resize(G->upperEdgeIdBound(), 0.0);

    G->parallelForEdges([&](node u, node v, edgeid eid) {
        Redundancy redundancy = getOverlap(u, v, neighbors, maxRank);

        scoreData[eid] = (double) redundancy.overlap;
    });

    hasRun = true;
}

} /* namespace NetworKit */
