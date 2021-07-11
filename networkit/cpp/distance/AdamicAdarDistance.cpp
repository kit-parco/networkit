// no-networkit-format
/*
 * AdamicAdarDistance.cpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/distance/AdamicAdarDistance.hpp>

namespace NetworKit {

AdamicAdarDistance::AdamicAdarDistance(const Graph& G) : NodeDistance(G) {
}

void AdamicAdarDistance::preprocess() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    Graph g = *G;

    //Node attribute: marker
    std::vector<bool> nodeMarker(g.upperNodeIdBound(), false);

    //Edge attribute: triangle count
    aaDistance = std::vector<double>(g.upperEdgeIdBound(), 0);

    g.forNodes([&](node u) {
        //Mark all neighbors
        g.forNeighborsOf(u, [&](node v) {
            nodeMarker[v] = true;
        });

        //For all neighbors: check for already marked neighbors.
        g.forNeighborsOf(u, [&](node, node v, edgeid eid_uv) {
            g.forNeighborsOf(v, [&](node, node w, edgeid eid_vw) {
                if (nodeMarker[w]) {

                    edgeid eid_uw = G->edgeId(u, w);

                    aaDistance[eid_uv] = aaDistance[eid_uv] + 1.0 / log(G->degree(w));
                    aaDistance[eid_uw] = aaDistance[eid_uw] + 1.0 / log(G->degree(v));
                    aaDistance[eid_vw] = aaDistance[eid_vw] + 1.0 / log(G->degree(u));
                }
            });

            nodeMarker[v] = false;
        });

        g.removeNode(u);
    });

    G->parallelForEdges([&](node, node, edgeid eid) {
        aaDistance[eid] = 1 / aaDistance[eid];
    });
}

double AdamicAdarDistance::distance(node u, node v) {
    edgeid eid = G->edgeId(u, v);
    return aaDistance[eid];
}


std::vector< double > AdamicAdarDistance::getEdgeScores() {
    return aaDistance;
}

} /* namespace NetworKit */
