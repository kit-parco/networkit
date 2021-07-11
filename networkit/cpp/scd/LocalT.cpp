#include <algorithm>
#include <limits>
#include <unordered_map>
#include <unordered_set>

#include <networkit/scd/LocalT.hpp>

#include "LocalDegreeDirectedGraph.hpp"

namespace NetworKit {

LocalT::LocalT(const Graph &g) : SelectiveCommunityDetector(g) {}

std::set<node> LocalT::expandOneCommunity(const std::set<node> &s) {
    // global data structures
    // result community
    std::set<node> result;

    // This algorithm creates a local graph. It contains the community and its shell.
    // stores the sum of the similarities of all nodes to nodes in the community
    std::vector<double> nodeInternalTriangles;
    std::vector<double> nodeExternalTriangles;
    std::vector<double> nodeSemiInternalTriangles;

    // indicates for every local node if it is in the community (true) or in the shell (false)
    std::vector<bool> inResult;

    std::vector<bool> inShell;

    count internalTriangles = 0, externalTriangles = 0;

    auto addNode = [&](node, node, double) {
        nodeInternalTriangles.push_back(0);
        nodeExternalTriangles.push_back(0);
        nodeSemiInternalTriangles.push_back(0);
        inResult.push_back(false);
        inShell.push_back(false);
    };

    LocalDegreeDirectedGraph<false, decltype(addNode)> localGraph(*g, addNode);

    std::unordered_set<node> shell;

    auto updateShell = [&](node u) {
        localGraph.forTrianglesOf(u, [&](node lv, node lw, edgeweight, edgeweight, edgeweight) {
            // collect counts and compute scores
            // new completely internal triangle, was not counted
            // previously
            if (inResult[lv] && inResult[lw]) {
                ++nodeInternalTriangles[lv];
                ++nodeInternalTriangles[lw];
                ++internalTriangles;
            } else if (inResult[lv] || inResult[lw]) {
                --externalTriangles;

                if (inResult[lv]) {
                    ++nodeInternalTriangles[lw];
                    --nodeSemiInternalTriangles[lw];
                } else {
                    ++nodeInternalTriangles[lv];
                    --nodeSemiInternalTriangles[lv];
                }
            } else {
                ++externalTriangles;
                if (inShell[lv]) {
                    --nodeExternalTriangles[lv];
                }
                ++nodeSemiInternalTriangles[lv];

                if (inShell[lw]) {
                    --nodeExternalTriangles[lw];
                }
                ++nodeSemiInternalTriangles[lw];
            }
        });

        g->forNeighborsOf(u, [&](node v) {
            node s = localGraph.ensureNodeExists(v);

            if (!inShell[s] && !inResult[s]) {
                shell.insert(s);
                inShell[s] = true;

                localGraph.forTrianglesOf(localGraph.toGlobal(s), [&](node lv, node lw, edgeweight,
                                                                      edgeweight, edgeweight) {
                    if (!inResult[lv] && !inResult[lw]) {
                        ++nodeExternalTriangles[s];
                    }
                });
            }
        });

#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
        count debugExternalTriangles = 0;
        count debugInternalTriangles = 0;
        for (node u : result) {
            count uExternal = 0, uInternal = 0;
            localGraph.forTrianglesOf(u, [&](node lv, node lw, edgeweight, edgeweight, edgeweight) {
                if (inResult[lv] && inResult[lw]) {
                    ++uInternal;
                } else if (!inResult[lv] && !inResult[lw]) {
                    ++uExternal;
                }
            });

            node lu = localGraph.ensureNodeExists(u);

            assert(uInternal == nodeInternalTriangles[lu]);
            assert(uInternal == nodeInternalTriangles[lu]);

            debugExternalTriangles += uExternal;
            debugInternalTriangles += uInternal;
        }

        assert(debugExternalTriangles == externalTriangles);
        assert(debugInternalTriangles / 3 == internalTriangles);

        for (node ls : shell) {
            count sExternal = 0, sInternal = 0, sSemiInternal = 0;
            localGraph.forTrianglesOf(localGraph.toGlobal(ls),
                                      [&](node lv, node lw, edgeweight, edgeweight, edgeweight) {
                                          if (inResult[lv] && inResult[lw]) {
                                              ++sInternal;
                                          } else if (inResult[lv] || inResult[lw]) {
                                              ++sSemiInternal;
                                          } else {
                                              ++sExternal;
                                          }
                                      });

            assert(sSemiInternal == nodeSemiInternalTriangles[ls]);
            assert(sInternal == nodeInternalTriangles[ls]);
            assert(sExternal == nodeExternalTriangles[ls]);
        }
#endif
#endif
    };

    // init community with seed set
    for (node u : s) {
        node lu = localGraph.ensureNodeExists(u);
        result.insert(u);
        shell.erase(lu);
        inResult[lu] = true;
        updateShell(u);
    }

    auto getScore = [](count intTriangles, count extTriangles) -> count {
        return std::max<int64_t>(
            0, intTriangles
                   * (static_cast<int64_t>(intTriangles) - static_cast<int64_t>(extTriangles)));
    };

    // expand (main loop)
    node uMax = none;
    do {

        uMax = none;
        count bestScore = getScore(internalTriangles, externalTriangles);
        count bestExternalTriangles = none;

        for (node lv : shell) {
            count newInternalTriangles = internalTriangles + nodeInternalTriangles[lv];
            count newExternalTriangles =
                externalTriangles + nodeExternalTriangles[lv] - nodeSemiInternalTriangles[lv];

            count newScore = getScore(newInternalTriangles, newExternalTriangles);

            if (newScore > bestScore
                || (newScore == bestScore && newExternalTriangles < bestExternalTriangles)) {
                uMax = lv;
                bestScore = newScore;
                bestExternalTriangles = newExternalTriangles;
            }
        }

        if (uMax != none) {
#ifdef NETWORKIT_SANITY_CHECKS
            assert(result.find(uMax) == result.end());
#endif
            node gUMax = localGraph.toGlobal(uMax);
            result.insert(gUMax);
            shell.erase(uMax);
            inResult[uMax] = true;
            updateShell(gUMax);
            assert(externalTriangles == bestExternalTriangles);
            assert(getScore(internalTriangles, externalTriangles) == bestScore);
        }

    } while (uMax != none);
    return result;
}

} /* namespace NetworKit */
