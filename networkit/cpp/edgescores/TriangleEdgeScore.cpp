// no-networkit-format
/*
 * TriangleEdgeScore.cpp
 *
 *  Created on: 29.08.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#include <omp.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/edgescores/TriangleEdgeScore.hpp>

namespace NetworKit {

TriangleEdgeScore::TriangleEdgeScore(const Graph& G) : EdgeScore<count>(G) {
}

void TriangleEdgeScore::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    // direct edge from high to low-degree nodes
    auto isOutEdge = [&](node u, node v) {
        return G->degree(u) > G->degree(v) || (G->degree(u) == G->degree(v) && u < v);
    };

    Aux::Timer filterEdgesTimer;
    filterEdgesTimer.start();
    // Store in-edges explicitly. Idea: all nodes have (relatively) low in-degree
    std::vector<index> inBegin(G->upperNodeIdBound() + 1);
    std::vector<node> inEdges(G->numberOfEdges());

    {
        index pos = 0;
        for (index u = 0; u < G->upperNodeIdBound(); ++u) {
            inBegin[u] = pos;
            if (G->hasNode(u)) {
                G->forEdgesOf(u, [&](node, node v, edgeid) {
                    if (isOutEdge(v, u)) {
                        inEdges[pos++] = v;
                    }
                });
            }
        }
        inBegin[G->upperNodeIdBound()] = pos;
    }

    filterEdgesTimer.stop();
    INFO("Needed ", filterEdgesTimer.elapsedMilliseconds(), "ms for filtering edges");

    //Edge attribute: triangle count
    std::vector<count> triangleCount(G->upperEdgeIdBound(), 0);
    // Store triangle counts of edges incident to the current node indexed by the adjacent node
    // none indicates that the edge to that node does not exist
    std::vector<std::vector<count> > incidentTriangleCount(omp_get_max_threads(), std::vector<count>(G->upperNodeIdBound(), none));

    Aux::Timer triangleTimer;
    triangleTimer.start();

    G->balancedParallelForNodes([&](node u) {
        auto tid = omp_get_thread_num();

        // mark nodes as neighbors
        G->forEdgesOf(u, [&](node, node v) {
            incidentTriangleCount[tid][v] = 0;
        });

        // Find all triangles of the form u-v-w-u where (v, w) is an in-edge.
        // Note that we find each triangle u is part of once.
        G->forEdgesOf(u, [&](node, node v) {
            // for all in-edges (v, w).
            for (index i = inBegin[v]; i < inBegin[v + 1]; ++i) {
                auto w = inEdges[i];

                if (incidentTriangleCount[tid][w] != none) {
                    // we have found a triangle u-v-w-u

                    // Record triangle count only for edges outgoing edges from u that should have consecutive edge ids
                    // This should be the same condition as in Graph::useEdgeInIteration and in Graph::indexEdges
                    // The count on (v, w) is updated when v or w is the central node

                    // Possibly record triangle for (u, v)
                    if (u >= v) {
                        ++incidentTriangleCount[tid][v];
                    }

                    // Possibly record triangle for (u, w)
                    if (u >= w) {
                        ++incidentTriangleCount[tid][w];
                    }
                }
            }
        });

        // Write local triangle counts into global array, unset local counters
        G->forEdgesOf(u, [&](node, node v, edgeid eid) {
            if (incidentTriangleCount[tid][v] > 0) {
                triangleCount[eid] += incidentTriangleCount[tid][v];
            }
            incidentTriangleCount[tid][v] = none;
        });
    });

    triangleTimer.stop();
    INFO("Needed ", triangleTimer.elapsedMilliseconds(), "ms for counting triangles");

    scoreData = std::move(triangleCount);
    hasRun = true;
}

count TriangleEdgeScore::score(node, node) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

count TriangleEdgeScore::score(edgeid) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
