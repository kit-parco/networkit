// no-networkit-format
/*
 * DynKatzCentrality.cpp
 *
 *  Created on: April 2018
 *      Author: Alexander van der Grinten
 *      based on code by Elisabetta Bergamini
 */

#include <cfloat>
#include <cmath>
#include <omp.h>

#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/centrality/DynKatzCentrality.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

DynKatzCentrality::DynKatzCentrality(const Graph& G, count k, bool groupOnly, double tolerance)
: Centrality(G, false), k(k), groupOnly(groupOnly), rankTolerance(tolerance),
        activeQueue(G.upperNodeIdBound()) {
    maxdeg = static_cast<count>(GraphTools::maxDegree(G));
    if (maxdeg == 0)
        throw std::runtime_error("Alpha is chosen based on the max. degree; therefore, that degree must not be zero");

    alpha = 1.0 / (maxdeg + 1.0);
    DEBUG("alpha: ", alpha);
    DEBUG("1/(1-alpha): ", 1/(1-alpha));
}

bool pairCompare(const std::pair<node, double>& firstElem, const std::pair<node, double>& secondElem) {
  return firstElem.second > secondElem.second;
}

void DynKatzCentrality::run() {
    DEBUG("DynKatz: Nodes: ", G.numberOfNodes(), ", edges: ", G.numberOfEdges(),
            ", graph is ", G.isDirected() ? "directed" : "undirected");

    nPaths.clear();
    nPaths.resize(1);
    nPaths[0].resize(G.upperNodeIdBound(), 1);

    isActive.clear();
    isActive.resize(G.upperNodeIdBound(), true);
    for(count u = 0; u < G.upperNodeIdBound(); u++)
        activeRanking.push_back(u);

    scoreData.clear();
    baseData.clear();
    boundData.clear();
    scoreData.resize(G.upperNodeIdBound(), 0.0);
    baseData.resize(G.upperNodeIdBound(), 0.0);
    boundData.resize(G.upperNodeIdBound(), DBL_MAX);

    levelReached = 0;
    do {
        doIteration();
    } while(!checkConvergence());
    DEBUG("DynKatz: Reached level: ", levelReached);

    hasRun = true;
}

void DynKatzCentrality::updateBatch(const std::vector<GraphEvent> &events){
    std::vector<count> preUpdatePaths(G.upperNodeIdBound(), 0);
    std::vector<count> preUpdateContrib(G.upperNodeIdBound(), 0);

    std::vector<bool> wasSeen(G.upperNodeIdBound(), false);
    std::vector<node> newlySeen;
    std::vector<node> seenNodes;

    count visitedEdges = 0;

    // First, we manually handle level 1. At level 1 only the two endpoints change.
    for(GraphEvent e : events) {
        if (e.type != GraphEvent::EDGE_ADDITION && e.type != GraphEvent::EDGE_REMOVAL) {
            throw std::runtime_error("Event type not allowed. Edge insertions or deletions only.");
        }

        if(!wasSeen[e.v]) {
            wasSeen[e.v] = true;
            seenNodes.push_back(e.v);
            preUpdatePaths[e.v] = 1;
        }

        if(!G.isDirected())
            if(!wasSeen[e.u]) {
                wasSeen[e.u] = true;
                seenNodes.push_back(e.u);
                preUpdatePaths[e.u] = 1;
            }
    }

    int maxThreads = omp_get_max_threads();

    count i = 1;
    while(i <= levelReached) {
        #pragma omp parallel for
        for (omp_index m = 0; m < static_cast<omp_index>(seenNodes.size()); ++m) {
            node v = seenNodes[m];
            preUpdateContrib[v] = preUpdatePaths[v];
            preUpdatePaths[v] = nPaths[i][v];
        }

        if(seenNodes.size() < G.numberOfNodes()) {
            // Caveat: We need preUpdateContrib[e.u] for deletions.
            // If e.u was not seen yet, we initialize that value here.
            for(GraphEvent e : events)
                if(e.type == GraphEvent::EDGE_REMOVAL && !wasSeen[e.u])
                    preUpdateContrib[e.u] = nPaths[i-1][e.u];

            // Subtract the old contribution and add the new one.
            for (node u : seenNodes) {
                // Note: For directed graphs here the direction has to be the opposite
                // of the static case.
                G.forEdgesOf(u, [&](node v, edgeweight) {
                    visitedEdges++;
                    if(!wasSeen[v]) {
                        wasSeen[v] = true;
                        newlySeen.push_back(v);
                        preUpdatePaths[v] = nPaths[i][v];
                    }

                    nPaths[i][v] -= preUpdateContrib[u];
                    nPaths[i][v] += nPaths[i-1][u];
                });
            }

            // If parallelism is available, it is faster to just recompute all nodes.
            // TODO: Check for parallelism.
            if(maxThreads > 1 && 2 * seenNodes.size() > G.numberOfNodes())
                G.forNodes([&] (node v) {
                    if(!wasSeen[v]) {
                        wasSeen[v] = true;
                        newlySeen.push_back(v);
                        preUpdatePaths[v] = nPaths[i][v];
                    }
                });

            // Handle the added/deleted edges.
            for(GraphEvent e : events) {
                if(e.type == GraphEvent::EDGE_ADDITION) {
                    nPaths[i][e.v] += nPaths[i-1][e.u];
                    if(!G.isDirected())
                        nPaths[i][e.u] += nPaths[i-1][e.v];
                }else{
                    assert(e.type == GraphEvent::EDGE_REMOVAL);
                    nPaths[i][e.v] -= preUpdateContrib[e.u];
                    if(!G.isDirected())
                        nPaths[i][e.u] -= preUpdateContrib[e.v];
                }
            }

            seenNodes.insert(seenNodes.end(), newlySeen.begin(), newlySeen.end());
            newlySeen.clear();

            // Update the Katz centrality from nPaths.
            auto alpha_pow = pow(alpha, i);
            #pragma omp parallel for
            for (omp_index m = 0; m < static_cast<omp_index>(seenNodes.size()); ++m) {
                node v = seenNodes[m];
                baseData[v] -= alpha_pow * preUpdatePaths[v];
                baseData[v] += alpha_pow * nPaths[i][v];
            }
        }else{
            // In this case, we're basically applying the static algorithm.
            auto alpha_pow = pow(alpha, i);
            G.balancedParallelForNodes([&](node u) {
                nPaths[i][u] = 0;
                G.forInEdgesOf(u, [&](node v, edgeweight) {
                    nPaths[i][u] += preUpdateContrib[v];
                });

                baseData[u] -= alpha_pow * preUpdatePaths[u];
                baseData[u] += alpha_pow * nPaths[i][u];
            });
        }

        i++;
    }

    i --;

    DEBUG("DynKatz: Done update iteration. visitedEdges = ", visitedEdges,
        ", speedup: ", double(levelReached*G.numberOfEdges()*2)/visitedEdges);

    // We compute the new bounds and reactive nodes here.

    // The following value is a lower bound that we know even without recomputing scoreData.
    double reactivation_threshold = 0.0;
    if(activeRanking.size() >= k)
        reactivation_threshold = baseData[*std::min_element(activeRanking.begin(),
                activeRanking.end(), [&] (node u, node v) {
            return baseData[u] < baseData[v];
        })];
    reactivation_threshold -= rankTolerance;
    DEBUG("DynKatz: Reactivation threshold: ", reactivation_threshold);

    auto alpha_pow = pow(alpha, i); // See doIteration().
    auto next_alpha_pow = alpha * alpha_pow;
    auto bound_factor = next_alpha_pow * (1/(1 - alpha * maxdeg));
    G.forNodes([&](node u){
        if(!G.isDirected()) {
            scoreData[u] = baseData[u] + next_alpha_pow * nPaths[i][u];
        }else{
            scoreData[u] = baseData[u];
        }
        boundData[u] = baseData[u] + nPaths[i][u] * bound_factor;

        // Reactivate nodes if they can potentially be in the top-k.
        if(!isActive[u] && boundData[u] >= reactivation_threshold) {
            isActive[u] = true;
            activeRanking.push_back(u);
        }
    });
    DEBUG("DynKatz: ", activeRanking.size(), " nodes in ranking after reactivation");

    // TODO what if the maxdeg increases???

    // Check if we need more iterations.
    if(!checkConvergence()) {
        do {
            doIteration();
        } while(!checkConvergence());
        DEBUG("DynKatz: Reached level: ", levelReached);
    }
}

double DynKatzCentrality::bound(node v) {
    assureFinished();
    return boundData.at(v);
}

bool DynKatzCentrality::areDistinguished(node u, node v) {
    if(scoreData[u] < scoreData[v])
        std::swap(u, v);
    return scoreData[u] > boundData[v];
}

bool DynKatzCentrality::areSufficientlyRanked(node high, node low) {
    return scoreData[high] > boundData[low] - rankTolerance;
}

void DynKatzCentrality::doIteration() {
    // The following variable is the level that his iteration will fill in.
    count r = levelReached + 1;

    nPaths.resize(r + 1);
    nPaths[r].resize(G.upperNodeIdBound(), 0);

    // Next, compute the ranking of active nodes for the current iteration.
    // GCC 6 is not smart enough to move the 'pow' out of the loop automatically.
    auto alpha_pow = pow(alpha, r);
    auto next_alpha_pow = alpha * alpha_pow;
    auto bound_factor = next_alpha_pow * (1/(1 - alpha * maxdeg));
    G.balancedParallelForNodes([&](node u){
        G.forInEdgesOf(u, [&](node v, edgeweight) {
            nPaths[r][u] += nPaths[r-1][v];
        });

        baseData[u] += alpha_pow * nPaths[r][u];
        // TODO: Enable this assertion.
//		assert(baseData[u] <= boundData[u]);
        if(!G.isDirected()) {
            scoreData[u] = baseData[u] + next_alpha_pow * nPaths[r][u];
        }else{
            scoreData[u] = baseData[u];
        }
        boundData[u] = baseData[u] + bound_factor * nPaths[r][u];
    });

    levelReached++;
}

bool DynKatzCentrality::checkConvergence() {
    if(useQueue) {
        if(!filledQueue) {
            G.forNodes([&] (node u) {
                activeQueue.insert(scoreData[u], u);
            });

            filledQueue = true;
        }

        // Rank the first node correctly.
        node v;
        do {
            v = activeQueue.peekMin(0).second;
            activeQueue.changeKey(scoreData[v], v);
        } while(v != activeQueue.peekMin(0).second);

        while(activeQueue.size() > 1) {
            // Rank the second node correctly.
            node u;
            do {
                u = activeQueue.peekMin(1).second;
                activeQueue.changeKey(scoreData[u], u);
            } while(u != activeQueue.peekMin(1).second);

            if(!areSufficientlyRanked(u, v))
                return false;
            activeQueue.extractMin();
            v = u;
        }

        return true;
    }else{
        // Deactivate nodes that cannot be in the top-k.
        if(activeRanking.size() > k) {
            // Doing only a partial sort here improves performance a lot.
            Aux::Timer sort_timer;
            sort_timer.start();
            std::partial_sort(activeRanking.begin(), activeRanking.begin() + k,
                    activeRanking.end(), [&] (node u, node v) {
                return scoreData[u] > scoreData[v];
            });
            sort_timer.stop();
            INFO("DynKatz: Partial sort time: ", sort_timer.elapsedMilliseconds());

            for(auto u : activeRanking)
                if(areSufficientlyRanked(activeRanking[k - 1], u))
                    isActive[u] = false;

            activeRanking.erase(std::remove_if(activeRanking.begin() + k,
                    activeRanking.end(), [&] (node u) {
                return !isActive[u];
            }), activeRanking.end());
        }

        assert(!activeRanking.empty());
/*
        double length = sqrt(G.parallelSumForNodes([&](node u) {
            return (scoreData[u] * scoreData[u]);
        }));
        DEBUG("DynKatz: Vector length: ", length);
*/
        DEBUG("DynKatz: In iteration ", levelReached, ": ", activeRanking.size(), " nodes remain.");

        if(activeRanking.size() > k)
            return false;

        // Once only k active nodes remain, we need to continue iterating
        // until the bounds seperate their position in the ranking.
        if(!groupOnly) {
            Aux::Timer sort_timer;
            sort_timer.start();
            std::sort(activeRanking.begin(), activeRanking.end(), [&] (node u, node v) {
                return scoreData[u] > scoreData[v];
            });
            sort_timer.stop();
            INFO("DynKatz: Sort time: ", sort_timer.elapsedMilliseconds());

            DEBUG("DynKatz: Node with highest centrality: ", activeRanking[0],
                    ", score: ", baseData[activeRanking[0]],
                    ", (upper) bound: ", boundData[activeRanking[0]],
                    ", lower bound: ", scoreData[activeRanking[0]]);

            for (count j = 1; j < std::min(size_t(k), activeRanking.size()); j ++) {
                node previous = activeRanking[j - 1];
                node current = activeRanking[j];
                if(!areSufficientlyRanked(previous, current))
                    return false;
            }
        }

        return true;
    }
}

} /* namespace NetworKit */
