// no-networkit-format
/*
 * ChibaNishizekiTriangleEdgeScore.cpp
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/edgescores/ChibaNishizekiTriangleEdgeScore.hpp>

namespace NetworKit {

ChibaNishizekiTriangleEdgeScore::ChibaNishizekiTriangleEdgeScore(const Graph& G) : EdgeScore<count>(G) {}

void ChibaNishizekiTriangleEdgeScore::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    std::vector<std::vector<std::pair<node, edgeid>>> edges(G->upperNodeIdBound());

    // copy edges with edge ids
    G->parallelForNodes([&](node u) {
        edges[u].reserve(G->degree(u));
        G->forEdgesOf(u, [&](node, node v, edgeid eid) {
            edges[u].emplace_back(v, eid);
        });
    });

    //Node attribute: marker
    std::vector<edgeid> nodeMarker(G->upperNodeIdBound(), none);

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
        //Mark all neighbors
        for (auto uv : edges[u]) {
            nodeMarker[uv.first] = uv.second;
        }

        //For all neighbors: check for already marked neighbors.
        for (auto uv : edges[u]) {
            auto& uedges = edges[uv.first];
            for (auto vw = uedges.begin(); vw != uedges.end(); ++vw) {
                // delete the edge to u as we do not need to consider it again.
                // the opposite edge doesn't need to be deleted as we will never again consider
                // outgoing edges of u as u cannot be reached anymore after the uv loop.
                if (vw->first == u) {
                    // move last element to current position in order to avoid changing too much
                    if (std::next(vw) == uedges.end()) {
                        uedges.pop_back();
                        break;
                    }

                    *vw = uedges.back();
                    uedges.pop_back();
                }

                if (nodeMarker[vw->first] != none) { // triangle found - count it!
                    edgeid eid_uw = nodeMarker[vw->first];

                    ++scoreData[uv.second];
                    ++scoreData[eid_uw];
                    ++scoreData[vw->second];
                }
            }

            nodeMarker[uv.first] = none; // all triangles with u and v have been counted already
        }
    }

    hasRun = true;
}

count ChibaNishizekiTriangleEdgeScore::score(node, node) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

count ChibaNishizekiTriangleEdgeScore::score(edgeid) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
