// no-networkit-format
/*
 * ChangeCorrectedTriangleScore.cpp
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#include <networkit/sparsification/ChanceCorrectedTriangleScore.hpp>

namespace NetworKit {

ChanceCorrectedTriangleScore::ChanceCorrectedTriangleScore(const Graph& G, const std::vector<count>& triangles) : EdgeScore<double>(G), triangles(&triangles) {
}

void ChanceCorrectedTriangleScore::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    scoreData.resize(G->upperEdgeIdBound(), 0.0);

    G->parallelForEdges([&](node u, node v, edgeid eid) {
        if ((*triangles)[eid] > 0) {
            scoreData[eid] = static_cast<double>((*triangles)[eid] * (G->numberOfNodes() - 2)) * 1.0 / static_cast<double>((G->degree(u) - 1) * (G->degree(v) - 1));
        } else if (G->degree(u) == 1 || G->degree(v) == 1) {
            scoreData[eid] = 1;
        }
    });

    hasRun = true;
}

double ChanceCorrectedTriangleScore::score(node, node) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

double ChanceCorrectedTriangleScore::score(edgeid) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
