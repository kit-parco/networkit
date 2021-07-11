// no-networkit-format
#include <networkit/centrality/LocalClusteringCoefficient.hpp>
#include <omp.h>

namespace NetworKit {

LocalClusteringCoefficient::LocalClusteringCoefficient(const Graph& G, bool turbo) : Centrality(G, false, false), turbo(turbo) {
    if (G.isDirected()) throw std::runtime_error("Not implemented: Local clustering coefficient is currently not implemented for directed graphs");
    if (G.numberOfSelfLoops()) throw std::runtime_error("Local Clustering Coefficient implementation does not support graphs with self-loops. Call Graph.removeSelfLoops() first.");
}

void LocalClusteringCoefficient::run() {
    count z = G.upperNodeIdBound();
    scoreData.clear();
    scoreData.resize(z); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

    std::vector<index> inBegin;
    std::vector<node> inEdges;

    if (turbo) {
        auto isOutEdge = [&](node u, node v) {
            return G.degree(u) > G.degree(v) || (G.degree(u) == G.degree(v) && u < v);
        };

        inBegin.resize(G.upperNodeIdBound() + 1);
        inEdges.resize(G.numberOfEdges());
        index pos = 0;
        for (index u = 0; u < G.upperNodeIdBound(); ++u) {
            inBegin[u] = pos;
            if (G.hasNode(u)) {
                G.forEdgesOf(u, [&](node v) {
                    if (isOutEdge(v, u)) {
                        inEdges[pos++] = v;
                    }
                });
            }
        }
        inBegin[G.upperNodeIdBound()] = pos;
    }

    std::vector<std::vector<bool> > nodeMarker(omp_get_max_threads());

    for (auto & nm : nodeMarker) {
        nm.resize(z, false);
    }


    G.balancedParallelForNodes([&](node u) {
        count d = G.degree(u);

        if (d < 2) {
            scoreData[u] = 0.0;
        } else {
            size_t tid = omp_get_thread_num();
            count triangles = 0;

            G.forEdgesOf(u, [&](node v) {
                nodeMarker[tid][v] = true;
            });

            G.forEdgesOf(u, [&](node, node v) {
                if (turbo) {
                    for (index i = inBegin[v]; i < inBegin[v+1]; ++i) {
                        node w = inEdges[i];
                        if (nodeMarker[tid][w]) {
                            triangles += 1;
                        }
                    }
                } else {
                    G.forEdgesOf(v, [&](node, node w) {
                        if (nodeMarker[tid][w]) {
                            triangles += 1;
                        }
                    });
                }
            });

            G.forEdgesOf(u, [&](node, node v) {
                nodeMarker[tid][v] = false;
            });

            scoreData[u] = (double) triangles / (double)(d * (d - 1)); // No division by 2 since triangles are counted twice as well!
            if (turbo) scoreData[u] *= 2; // in turbo mode, we count each triangle only once
        }
    });
    hasRun = true;
}


double LocalClusteringCoefficient::maximum() {
    return 1.0;
}

}
