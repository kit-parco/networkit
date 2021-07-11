// no-networkit-format
/*
 * ClusteringCoefficient.cpp
 *
 *  Created on: 08.04.2013
 *      Author: Lukas Barth, David Weiss
 */

#include <omp.h>
#include <unordered_set>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/centrality/LocalClusteringCoefficient.hpp>
#include <networkit/global/ClusteringCoefficient.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

double ClusteringCoefficient::sequentialAvgLocal(const Graph &G) {
    WARN("DEPRECATED: use centrality.LocalClusteringCoefficient and take average");
    std::vector<std::vector<node> > edges(G.upperNodeIdBound());

    // copy edges with edge ids
    G.parallelForNodes([&](node u) {
        edges[u].reserve(G.degree(u));
        G.forEdgesOf(u, [&](node, node v, edgeid) {
            edges[u].emplace_back(v);
        });
    });

    //Node attribute: marker
    std::vector<bool> nodeMarker(G.upperNodeIdBound(), false);

    //Edge attribute: triangle count
    std::vector<count> triangleCount(G.upperNodeIdBound(), 0);

    // bucket sort
    count n = G.numberOfNodes();
    std::vector<node> sortedNodes(n);
    {
        std::vector<index> nodePos(n + 1, 0);

        G.forNodes([&](node u) {
            ++nodePos[n - G.degree(u)];
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

        G.forNodes([&](node u) {
            sortedNodes[nodePos[n - G.degree(u)]++] = u;
        });
    }

    for (node u : sortedNodes) {
        //Mark all neighbors
        for (auto v : edges[u]) {
            nodeMarker[v] = true;
        }

        //For all neighbors: check for already marked neighbors.
        for (auto v : edges[u]) {
            for (auto w = edges[v].begin(); w != edges[v].end(); ++w) {
                // delete the edge to u as we do not need to consider it again.
                // the opposite edge doesn't need to be deleted as we will never again consider
                // outgoing edges of u as u cannot be reached anymore after the uv loop.
                if (*w == u) {
                    // move last element to current position in order to avoid changing too much
                    *w = edges[v].back();
                    edges[v].pop_back();
                    if (w == edges[v].end()) // break if we were at the last element already
                        break;
                }

                if (nodeMarker[*w]) { // triangle found - count it!
                    ++triangleCount[u];
                    ++triangleCount[*w];
                    ++triangleCount[v];
                }
            }

            nodeMarker[v] = false; // all triangles with u and v have been counted already
        }
    }

    double coefficient = 0;
    count size = 0;
    G.forNodes([&](node u) {
        count d = G.degree(u);
        if (d > 1) {
            coefficient += static_cast<double>(triangleCount[u]) * 2.0 / static_cast<double>(d * (d - 1));
            size++;
        }
    });

    if (size == 0) {
        return 0; // no triangle exists
    }

    return coefficient / size;
}

double ClusteringCoefficient::avgLocal(Graph& G, bool turbo) {
    WARN("DEPRECATED: use centrality.LocalClusteringCoefficient and take average");
    LocalClusteringCoefficient lcc(G, turbo);
    lcc.run();
    auto coefficients = lcc.scores(); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

    double sum = 0.0;
    count size = 0;

    G.forNodes([&](node u) {
        if (G.degree(u) >= 2) {
            sum += coefficients[u];
            size++;
        }
    });

    if (size == 0) {
        return 0; // no triangle exists
    }

    return sum / (double) size;
}

double ClusteringCoefficient::approxAvgLocal(Graph& G, const count trials) {

    double triangles = 0;
    for (count k = 0; k < trials; ++k) {
        node v = GraphTools::randomNode(G);
        TRACE("trial ", k, " sampled node ", v);

        if (G.degree(v) < 2) {
            // this vertex can never be part of a triangle,
            // nor middle point of a path of length 3
            --k;  // do not count trial
            continue;
        }

        TRACE("deg(v) = ", G.degree(v));
        node u = GraphTools::randomNeighbor(G, v);
        node w = GraphTools::randomNeighbor(G, v);
        TRACE("u=", u);
        TRACE("w=", w);

        // TODO This could be sped up for degree(v) == 2...
        while (u == w) {
            w = GraphTools::randomNeighbor(G, v);
            TRACE("w=", w);
        }

        if (G.hasEdge(u,w)) {
            triangles++;
        }
    }

    return triangles / (double) trials;
}


double ClusteringCoefficient::exactGlobal(Graph& G) {
    count z = G.upperNodeIdBound();
    std::vector<count> triangles(z); // triangles including node u (every triangle is counted six times)

    std::vector<std::vector<bool> > nodeMarker(omp_get_max_threads());
    for (auto &nm : nodeMarker) {
        nm.resize(z, false);
    }

    G.balancedParallelForNodes([&](node u){

        size_t tid = omp_get_thread_num();
        count tr = 0;

        if (G.degree(u) > 1) {
            G.forEdgesOf(u, [&](node, node v) {
                nodeMarker[tid][v] = true;
            });

            G.forEdgesOf(u, [&](node, node v) {
                G.forEdgesOf(v, [&](node, node w) {
                    if (nodeMarker[tid][w]) {
                        tr += 1;
                    }
                });
            });

            G.forEdgesOf(u, [&](node, node v) {
                nodeMarker[tid][v] = false;
            });
        }

        triangles[u] = tr;
    });

  double denominator = G.parallelSumForNodes([&](node u){
        return G.degree(u) * (G.degree(u) - 1);
    });

    double cc = G.parallelSumForNodes([&](node u){
        return triangles[u];
    });

    if (denominator == 0) {
        return 0; // no triangle exists
    }

    cc /= denominator;

    return cc;
}


double ClusteringCoefficient::approxGlobal(Graph& G, const count trials) {
    count z = G.upperNodeIdBound();

  // Calculate prefix sum over the nodes where each node v counts deg(v)*(deg(v)-1) times
    std::vector<count> weight(z);
    count psum = 0;
    G.forNodes([&](node v) {
        psum += G.degree(v) * (G.degree(v) - 1);
        weight[v] = psum;
    });

    if (psum == 0) return 0; // no node has degree > 1 - no triangle exists!

    // WARNING: I assume RAND_MAX to be larger than PSUM. If this should not hold for an application
    // or implementation of the standard library, a more sophisticated version of determining a
    // vertex uniformly at random must be used.

    double triangles = 0;
    for (count k = 0; k < trials; ++k) {
        count r = Aux::Random::integer(psum - 1);

        // plain old binary search:
        index low = 0;
        index high = G.upperNodeIdBound();
        while (low < high) {
            index middle = (low + high) / 2;

            if (weight[middle] < r) {
                low = middle + 1;
            } else if (weight[middle] > r) {
                high = middle;
            } else {
                low = high = middle;
            }
        }

        node v = low; // ewww.. setting a vertex to an index.. but works.

        count vDeg = G.degree(v);
        if (vDeg < 2) {
            // this vertex can never be part of a triangle,
            // nor middle point of a path of length 3
            --k; // do not count trial
            continue;
        }

        node u = GraphTools::randomNeighbor(G, v);
        node w = GraphTools::randomNeighbor(G, v);

        // TODO This could be sped up for degree(v) == 2...
        while (u == w) {
            w = GraphTools::randomNeighbor(G, v);
        }

        if (G.hasEdge(u,w)) {
            triangles++;
        }
    }

    return triangles / trials;
}

} /* namespace NetworKit */
