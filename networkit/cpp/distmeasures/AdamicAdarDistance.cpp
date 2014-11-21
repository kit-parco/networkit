/*
 * AdamicAdarDistance.cpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#include "AdamicAdarDistance.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Timer.h"

namespace NetworKit {

AdamicAdarDistance::AdamicAdarDistance(const Graph& G) : NodeDistance(G) {
}

void AdamicAdarDistance::preprocess() {
	Graph g = G;

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
		g.forNeighborsOf(u, [&](node _u, node v, edgeid eid_uv) {
			g.forNeighborsOf(v, [&](node _v, node w, edgeid eid_vw) {
				if (nodeMarker[w]) {

					edgeid eid_uw = G.edgeId(u, w);

					aaDistance[eid_uv] = aaDistance[eid_uv] + 1.0 / log(g.degree(w));
					aaDistance[eid_uw] = aaDistance[eid_uw] + 1.0 / log(g.degree(v));
					aaDistance[eid_vw] = aaDistance[eid_vw] + 1.0 / log(g.degree(u));
				}
			});

			nodeMarker[v] = false;
		});

		removeNode(g, u);
	});

	G.parallelForEdges([&](node u, node v, edgeid eid) {
		aaDistance[eid] = 1 / aaDistance[eid];
	});
}

double AdamicAdarDistance::distance(node u, node v) {
	edgeid eid = G.edgeId(u, v);
	return aaDistance[eid];
}


std::vector< double > AdamicAdarDistance::getEdgeAttribute() {
	return aaDistance;
}

void AdamicAdarDistance::removeNode(Graph& graph, node u) {
	//isolate the node before removing it.
	graph.forNeighborsOf(u, [&](node v) {
		graph.removeEdge(u,v);
	});

	graph.removeNode(u);
}

} /* namespace NetworKit */
