// no-networkit-format
/*
 * JaccardDistance.cpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/distance/JaccardDistance.hpp>

namespace NetworKit {

JaccardDistance::JaccardDistance(const Graph& G, const std::vector<count>& triangles) : NodeDistance(G), triangles(triangles) {
}

void JaccardDistance::preprocess() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    jDistance = std::vector<double>(G->upperEdgeIdBound());

    G->parallelForEdges([&](node u, node v, edgeid eid) {
        jDistance[eid] = getJaccardDistance(G->degree(u), G->degree(v), triangles[eid]);
    });
}

double JaccardDistance::distance(node u, node v) {
    edgeid eid = G->edgeId(u, v);
    return getJaccardDistance(G->degree(u), G->degree(v), triangles[eid]);
}

std::vector<double> JaccardDistance::getEdgeScores() {
    return jDistance;
}

inline double JaccardDistance::getJaccardDistance(count degU, count degV, count t) {
    return 1.0 - (static_cast<double>(t * 1.0) / static_cast<double>(degU + degV - t));
}

} /* namespace NetworKit */
