#include <algorithm>
#include <limits>
#include <unordered_map>

#include <tlx/container/d_ary_addressable_int_heap.hpp>
#include <tlx/unused.hpp>

#include <networkit/scd/TCE.hpp>

#include "LocalDegreeDirectedGraph.hpp"

namespace NetworKit {

TCE::TCE(const Graph &g, bool refine, bool useJaccard)
    : SelectiveCommunityDetector(g), refine(refine), useJaccard(useJaccard) {}

namespace {

#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
/*
 * Only used for debugging.
 */
double weightedEdgeScore(const Graph &g, node, node v, edgeweight uvWeight, double uDegree,
                         const std::unordered_map<node, double> &uNeighbors, bool useJaccard) {
    double nom = uvWeight;
    double vDegree = 0.0;

    g.forNeighborsOf(v, [&](node, node w, edgeweight vwWeight) {
        vDegree += vwWeight;
        auto uw = uNeighbors.find(w);
        if (uw != uNeighbors.end())
            nom += std::min(vwWeight, uw->second);
    });

    if (vDegree == 0.0)
        return 0.0;

    double denom = useJaccard ? (uDegree + vDegree - nom) : std::min(uDegree, vDegree);
    return nom / (denom * g.degree(v));
}
#endif
#endif

template <bool isWeighted>
std::set<node> expandSeedSetInternal(const Graph &g, const std::set<node> &s, bool refine,
                                     bool useJaccard) {
    // global data structures
    std::set<node> result = s;

    std::vector<double> weightedDegrees;
    std::vector<double> nodeScore;

    std::vector<double> cutEdges;
    std::vector<bool> inResult;

    // heap that contains the nodes of the shell that still need to be considered
    class Compare {
        const std::vector<double> &similarity;

    public:
        Compare(const std::vector<double> &similarity) : similarity(similarity) {}

        bool operator()(node u, node v) const { return similarity[u] > similarity[v]; }
    };

    tlx::d_ary_addressable_int_heap<node, 4, Compare> shell((Compare(nodeScore)));

    // data structures only for neighbors of a node
    std::vector<double> triangleSum;

    auto nodeAdded = [&](node, node, double weightedDegree) {
        weightedDegrees.push_back(weightedDegree);
        nodeScore.push_back(0);
        cutEdges.push_back(0);

        inResult.push_back(false);
        triangleSum.push_back(0);
    };

    LocalDegreeDirectedGraph<isWeighted, decltype(nodeAdded)> localGraph(g, nodeAdded);

    auto updateShell = [&](node u, node lu, bool noverify) -> double {
#if defined(NDEBUG) || !defined(NETWORKIT_SANITY_CHECKS)
        tlx::unused(noverify); // only used in debug mode
#endif
        if (g.degree(u) == 0)
            return 0.0;

        double xDegree = weightedDegrees[lu];

        localGraph.forTrianglesOf(
            u, [&](node lv, node y, double weightUv, double weightUy, double weightVy) {
                if (isWeighted) {
                    triangleSum[y] += std::min(weightUv, weightVy);
                    triangleSum[lv] += std::min(weightUy, weightVy);
                } else {
                    triangleSum[y] += 1;
                    triangleSum[lv] += 1;
                }
            });

#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
        std::unordered_map<node, double> uNeighbors;
        g.forNeighborsOf(
            u, [&](node, node v, edgeweight ew) { uNeighbors.insert(std::make_pair(v, ew)); });
#endif
#endif

        // collect counts and compute scores
        localGraph.forLocalNeighbors([&](node v, edgeweight weight) {
            if (!inResult[v]) {
                double nom = weight + triangleSum[v];
                double wd = weightedDegrees[v];

                double scoreUv = 0.0;
                if (wd > 0.0) {
                    double denom = useJaccard ? (wd + xDegree - nom) : std::min(wd, xDegree);
                    scoreUv = nom / (denom * g.degree(localGraph.toGlobal(v)));
                }

#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
                assert(scoreUv
                       == weightedEdgeScore(g, u, localGraph.toGlobal(v), weight, xDegree,
                                            uNeighbors, useJaccard));
#endif
#endif
                nodeScore[v] += scoreUv;
                shell.update(v);

                cutEdges[v] += weight;

#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
                if (!noverify) {
                    double debugWeight = 0;
                    g.forNeighborsOf(localGraph.toGlobal(v), [&](node, node y, edgeweight ew) {
                        if (result.count(y)) {
                            debugWeight += ew;
                        }
                    });
                    assert(debugWeight == cutEdges[v]);
                }
#endif
#endif
            }

            // reset triangle_sum
            triangleSum[v] = 0;
        });

#ifdef NETWORKIT_SANITY_CHECKS
        assert(
            std::all_of(triangleSum.begin(), triangleSum.end(), [](double s) { return s == 0; }));
#endif

        return xDegree;
    };

    // init community with seed set
    double volume = 0.0;
    double numCutEdges = 0.0;

    for (auto u : result) {
        node lu = localGraph.ensureNodeExists(u);
        inResult[lu] = true;
    }

    for (auto u : result) {
        volume += updateShell(u, localGraph.ensureNodeExists(u), true);
    }

    for (const auto &vv : cutEdges) {
        numCutEdges += vv;
    }

    // expand (main loop)
    while (!shell.empty()) {
        node uMax = shell.extract_top();
        node gUMax = localGraph.toGlobal(uMax);
        double uMaxVol = weightedDegrees[uMax];

        double numCutEdgesNew = numCutEdges + uMaxVol - (2 * cutEdges[uMax]);
        double volumeNew = volume + uMaxVol;

#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
        {
            double debugVolumeOld = 0.0;
            double debugCutsizeOld = 0.0;

            double debugVolumeNew = 0.0;
            double debugCutsizeNew = 0.0;

            for (auto u : result) {
                g.forEdgesOf(u, [&](node, node v, edgeweight ew) {
                    debugVolumeOld += ew;
                    debugVolumeNew += ew;
                    if (!result.count(v)) {
                        debugCutsizeOld += ew;

                        if (v != gUMax) {
                            debugCutsizeNew += ew;
                        }
                    }
                });
            }
            g.forEdgesOf(gUMax, [&](node, node v, edgeweight ew) {
                debugVolumeNew += ew;
                if (!result.count(v)) {
                    debugCutsizeNew += ew;
                }
            });

            assert(debugVolumeOld == volume);
            assert(debugVolumeNew == volumeNew);
            assert(debugCutsizeNew == numCutEdgesNew);
            assert(debugCutsizeOld == numCutEdges);
        }
#endif
#endif

        // if conductance decreases (i.e., improves)
        if ((numCutEdgesNew / volumeNew) < (numCutEdges / volume)) {
            result.insert(gUMax);
            inResult[uMax] = true;
            updateShell(gUMax, uMax, false);

            numCutEdges = numCutEdgesNew;
            volume = volumeNew;
        }
    }

    if (refine) {
        for (auto it = result.begin(); it != result.end();) {
            node u = *it;
            double uVol = 0;
            double uCutChange = 0;

            g.forNeighborsOf(u, [&](node, node v, edgeweight ew) {
                uVol += ew;
                if (result.count(v) > 0) {
                    uCutChange += ew;
                } else {
                    uCutChange -= ew;
                }
            });

            double numCutEdgesNew = numCutEdges + uCutChange;
            double volumeNew = volume - uVol;

#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
            {
                double debugVolumeOld = 0.0;
                double debugCutsizeOld = 0.0;

                double debugVolumeNew = 0.0;
                double debugCutsizeNew = 0.0;

                for (auto v : result) {
                    g.forEdgesOf(v, [&](node, node x, edgeweight ew) {
                        debugVolumeOld += ew;
                        if (v != u) {
                            debugVolumeNew += ew;
                        }
                        if (!result.count(x)) {
                            debugCutsizeOld += ew;
                            if (v != u) {
                                debugCutsizeNew += ew;
                            }
                        }
                        if (x == u) {
                            debugCutsizeNew += ew;
                        }
                    });
                }

                assert(debugVolumeOld == volume);
                assert(debugVolumeNew == volumeNew);
                assert(debugCutsizeNew == numCutEdgesNew);
                assert(debugCutsizeOld == numCutEdges);
            }
#endif
#endif

            if ((numCutEdgesNew / volumeNew) < (numCutEdges / volume)) {
                numCutEdges = numCutEdgesNew;
                volume = volumeNew;

                TRACE("Removing ", u, " again");

                it = result.erase(it);
            } else {
                ++it;
            }
        }
    }

    return result;
}

} // namespace

std::set<node> TCE::expandOneCommunity(const std::set<node> &s) {
    if (g->isWeighted()) {
        return expandSeedSetInternal<true>(*g, s, refine, useJaccard);
    } else {
        return expandSeedSetInternal<false>(*g, s, refine, useJaccard);
    }
}

} /* namespace NetworKit */
