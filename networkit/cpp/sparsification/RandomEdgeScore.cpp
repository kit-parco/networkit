// no-networkit-format
/*
 * RandomEdgeScore.cpp
 *
 *  Created on: 11.08.2014
 *      Author: Gerd Lindner
 */

#include <networkit/sparsification/RandomEdgeScore.hpp>

namespace NetworKit {

RandomEdgeScore::RandomEdgeScore(const Graph& G) : EdgeScore<double>(G) {}

void RandomEdgeScore::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }
    scoreData.resize(G->upperEdgeIdBound(), 0.0);

    G->parallelForEdges([&](node, node, edgeid eid) {
        scoreData[eid] = Aux::Random::probability();
    });
    hasRun = true;
}

double RandomEdgeScore::score(node, node) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

double RandomEdgeScore::score(edgeid) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
