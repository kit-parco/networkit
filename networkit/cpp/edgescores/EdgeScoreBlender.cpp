// no-networkit-format
/*
 * EdgeScoreBlender.cpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#include <networkit/edgescores/EdgeScoreBlender.hpp>

namespace NetworKit {

EdgeScoreBlender::EdgeScoreBlender(const Graph &G, const std::vector< double > &attribute0, const std::vector< double > &attribute1, const std::vector< bool > &selection) :
EdgeScore(G), attribute0(&attribute0), attribute1(&attribute1), selection(&selection) {}

void EdgeScoreBlender::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    scoreData.resize(G->upperEdgeIdBound());

    G->parallelForEdges([&](node, node, edgeid eid) {
        scoreData[eid] = ((*selection)[eid] ? (*attribute1)[eid] : (*attribute0)[eid]);
    });

    hasRun = true;
}

double EdgeScoreBlender::score(node, node) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

double EdgeScoreBlender::score(edgeid) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
