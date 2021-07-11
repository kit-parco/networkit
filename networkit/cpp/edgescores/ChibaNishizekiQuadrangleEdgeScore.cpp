// no-networkit-format
/*
 * ChibaNishizekiQuadrangleEdgeScore.cpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#include <networkit/edgescores/ChibaNishizekiQuadrangleEdgeScore.hpp>

namespace NetworKit {

ChibaNishizekiQuadrangleEdgeScore::ChibaNishizekiQuadrangleEdgeScore(const Graph& G) : EdgeScore<count>(G) {}

void ChibaNishizekiQuadrangleEdgeScore::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    std::vector<std::vector<std::pair<node, edgeid> > > edges(G->upperNodeIdBound());

    // copy edges with edge ids
    G->parallelForNodes([&](node u) {
        edges[u].reserve(G->degree(u));
        G->forEdgesOf(u, [&](node, node v, edgeid eid) {
            edges[u].emplace_back(v, eid);
        });
    });

    //Node attribute: marker
    std::vector<count> nodeMarker(G->upperNodeIdBound(), 0);

    //Edge attribute: triangle count
    scoreData.resize(G->upperEdgeIdBound(), 0);

    // bucket sort
    count n = G->numberOfNodes();
    std::vector<node> sortedNodes(n);
    {
        std::vector<index> nodePos(n + 1, 0);

        G->forNodes([&](node u) {
            ++nodePos[n - G->degree(u)];
        });

        // exclusive prefix sum
        index tmp = nodePos[0];
        index sum = tmp;
        nodePos[0] = 0;

        for (index i = 1; i < nodePos.size(); ++i) {
            tmp = nodePos[i];
            nodePos[i] = sum;
            sum += tmp;
        }

        G->forNodes([&](node u) {
            sortedNodes[nodePos[n - G->degree(u)]++] = u;
        });
    }

    for (node u : sortedNodes) {
        for (auto uv : edges[u]) {
            for (auto vw = edges[uv.first].begin(); vw != edges[uv.first].end(); ++vw) {
                // delete the edge to u as we do not need to consider it again.
                // the opposite edge doesn't need to be deleted as we will never again consider
                // outgoing edges of u as u cannot be reached anymore after the uv loop.
                if (vw->first == u) {
                    if (std::next(vw) == edges[uv.first].end()) {
                        edges[uv.first].pop_back();
                        break;
                    }

                    // move last element to current position in order to avoid changing too much
                    *vw = edges[uv.first].back();
                    edges[uv.first].pop_back();
                }

                ++nodeMarker[vw->first];
            }
        }

        for (auto uv : edges[u]) {
            for (auto vw : edges[uv.first]) {
                if (nodeMarker[vw.first] > 1) {
                    scoreData[uv.second] += nodeMarker[vw.first] - 1;
                    scoreData[vw.second] += nodeMarker[vw.first] - 1;
                }
            }
        }

        for (auto uv : edges[u]) {
            for (auto vw : edges[uv.first]) {
                nodeMarker[vw.first] = 0;
            }
        }
    }

    hasRun = true;
}

count ChibaNishizekiQuadrangleEdgeScore::score(node, node) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

count ChibaNishizekiQuadrangleEdgeScore::score(edgeid) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

}/* namespace NetworKit */
