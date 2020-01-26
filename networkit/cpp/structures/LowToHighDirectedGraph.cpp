/*
 * LowToHighDirectedGraph.cpp
 *
 * Created: 2018-12-12
 * Author: Armin Wiebigke
 */

#include "networkit/structures/LowToHighDirectedGraph.hpp"

namespace NetworKit {

LowToHighDirectedGraph::LowToHighDirectedGraph(const NetworKit::Graph &G) {
    edgesBegin.resize(G.upperNodeIdBound() + 1);
    edgeTargets.resize(G.numberOfEdges());
    if (G.isWeighted()) {
        edgeWeights.resize(G.numberOfEdges());
    }

    // direct edge from low to high-degree nodes
    auto isOutEdge = [&](node u, node v) {
        return G.degree(u) < G.degree(v) || (G.degree(u) == G.degree(v) && u < v);
    };

    index pos = 0;
    for (node u = 0; u < G.upperNodeIdBound(); ++u) {
        edgesBegin[u] = pos;
        if (G.hasNode(u)) {
            G.forEdgesOf(u, [&](node, node v, edgeweight ew) {
                if (isOutEdge(u, v)) {
                    if (G.isWeighted()) {
                        edgeWeights[pos] = ew;
                    }
                    edgeTargets[pos++] = v;
                }
            });
        }
    }

    edgesBegin[G.upperNodeIdBound()] = pos;
}

std::vector<node> LowToHighDirectedGraph::neighborsOf(node u) const {
    std::vector<node> neighbors(edgesBegin[u + 1] - edgesBegin[u]);
    std::copy(edgeTargets.begin() + edgesBegin[u], edgeTargets.begin() + edgesBegin[u + 1],
              neighbors.begin());
    return neighbors;
}

} /* namespace NetworKit */
