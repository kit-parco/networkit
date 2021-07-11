// no-networkit-format
/*
 * TopCloseness.cpp
 *
 *  Created on: 28.02.2018
 *      Author: nemes, Eugenio Angriman
 */

#include <cmath>
#include <limits>
#include <omp.h>

#include <networkit/auxiliary/PrioQueue.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/distance/AffectedNodes.hpp>
#include <networkit/centrality/DynTopHarmonicCloseness.hpp>
#include <networkit/graph/BFS.hpp>

namespace NetworKit {

DynTopHarmonicCloseness::DynTopHarmonicCloseness(const Graph &G, count k,
                                                 bool useBFSbound)
    : G(G), k(k), useBFSbound(useBFSbound),
      allScores(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
      isExact(G.upperNodeIdBound(), false),
      isValid(G.upperNodeIdBound(), false),
      cutOff(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
      exactCutOff(G.upperNodeIdBound(), 0), hasComps(false),
      component(G.upperNodeIdBound()), rOld(G.upperNodeIdBound()) {}

DynTopHarmonicCloseness::~DynTopHarmonicCloseness() {
    if (hasComps) {
        delete comps;
    }
    if (hasWComps) {
        delete wComps;
    }
}

std::pair<edgeweight, bool>
DynTopHarmonicCloseness::BFScut(node v, edgeweight x, count n, count r,
                                std::vector<uint8_t> &visited,
                                std::vector<count> &distances,
                                std::vector<node> &pred, count &visitedEdges) {

    double d = 0;
    int64_t gamma = 0, nd = 0;

    distances[v] = 0;
    std::queue<node> Q;
    std::queue<node> toReset;
    Q.push(v);
    toReset.push(v);
    ++nd;

    visited[v] = true;

    /* Clean up the visited vector for the next run.
     * We don't need to do this for the distances and pred vector
     * because those are only read after visited[u] is true.
     */
    auto cleanup = [&]() {
        while (!toReset.empty()) {
            node w = toReset.front();
            toReset.pop();
            visited[w] = false;
        }
    };

    edgeweight c = 0;
    edgeweight ctilde = edgeweight(n - 1);

    do {
        node u = Q.front();
        Q.pop();

        if (distances[u] > d) {
            d = d + 1;

            double d2 = static_cast<double>(int64_t(r) - nd);
            ctilde = c + (static_cast<edgeweight>(gamma) / ((d + 2.0) * (d + 1.0))) +
                     (d2 / (d + 2.0));

            if (ctilde < x) {
                exactCutOff[v] = true;
                cutOff[v] = d;
                cleanup();
                return std::make_pair(ctilde, false);
            }

            gamma = 0;
        }
        bool cont = true;
        G.forNeighborsOf(u, [&](node w) {
            if (cont) {
                visitedEdges++;
                if (!visited[w]) {
                    distances[w] = distances[u] + 1;

                    Q.push(w);
                    toReset.push(w);
                    visited[w] = true;

                    c += 1.0 / distances[w];

                    if (!G.isDirected()) {
                        gamma += (G.degree(w) - 1); // notice: only because it's undirected
                    } else {
                        gamma += G.degree(w);
                    }

                    ++nd;
                    pred[w] = u;
                } else if (distances[w] > 1 && (pred[u] != w)) {
                    ctilde = ctilde - 1.0 / (d + 1.0) + 1.0 / (d + 2.0);

                    if (ctilde < x) {
                        cont = false;
                    }
                }
            }
        });
        if (ctilde < x) {
            exactCutOff[v] = false;
            cutOff[v] = d;
            cleanup();

            return std::make_pair(ctilde, false);
        }
    } while (!Q.empty());

    cleanup();
    exactCutOff[v] = false;
    cutOff[v] = std::numeric_limits<edgeweight>::max();

    return std::make_pair(c, true);
}

void DynTopHarmonicCloseness::BFSbound(node source, std::vector<double> &S2,
                                       count *visEdges) {
    count r = 0;
    n = G.upperNodeIdBound();
    std::vector<std::vector<node>> levels(n);
    std::vector<count> nodesPerLev(n, 0);
    std::vector<count> sumLevs(n, 0);
    count nLevs = 0;
    levels[nLevs].clear();
    double sum_dist = 0;
    edgeweight level_bound;

    auto inverseDistance = [&](edgeweight dist) { return 1.0 / dist; };

    Traversal::BFSfrom(G, source, [&](node u, count dist) {
        sum_dist += dist > 0 ? inverseDistance(dist) : 0;

        ++r;
        if (dist > nLevs) {
            sumLevs[nLevs] += nodesPerLev[nLevs];
            sumLevs[nLevs + 1] = sumLevs[nLevs];
            ++nLevs;
            levels[nLevs].clear();
        }
        levels[nLevs].push_back(u);
        nodesPerLev[nLevs]++;
    });
    sumLevs[nLevs] += nodesPerLev[nLevs];
    if (G.isDirected()) {
        (*visEdges) += G.numberOfEdges();
    } else {
        (*visEdges) += 2 * G.numberOfEdges();
    }

    S2[source] = sum_dist;

    // we compute the bound for the first level
    double closeNodes = 0, farNodes = 0;
    for (count j = 0; j <= nLevs; j++) {
        if (j <= 2) {
            closeNodes += nodesPerLev[j];
        } else {
            farNodes +=
                nodesPerLev[j] * inverseDistance(double(std::abs((double)j - 1.)));
        }
    }

    level_bound = inverseDistance(2.) * closeNodes + farNodes;

    for (count j = 0; j < levels[1].size(); j++) {
        node w = levels[1][j];
        // we subtract 2 not to count the node itself
        double bound = (level_bound - inverseDistance(2.) +
                        (inverseDistance(1.) - inverseDistance(2.)) * G.degree(w));

        if (bound < S2[w] &&
            (!G.isDirected() || component[w] == component[source])) {
            S2[w] = bound;
        }
    }

    // now we compute it for the other levels
    for (count i = 2; i <= nLevs; i++) {

        level_bound = 0;
        // TODO: OPTIMIZE?
        if (!G.isDirected()) {
            for (count j = 0; j <= nLevs; j++) {
                level_bound += inverseDistance(std::max(
                                   2., double(std::abs((double)j - (double)i)))) *
                               nodesPerLev[j];
            }
        } else {
            for (count j = 0; j <= nLevs; j++) {
                level_bound += inverseDistance(std::max(2., (double)j - (double)i)) *
                               nodesPerLev[j];
            }
        }

        for (count j = 0; j < levels[i].size(); ++j) {
            node w = levels[i][j];
            double bound =
                (level_bound - inverseDistance(2.) +
                 (inverseDistance(1.) - inverseDistance(2.)) * G.degree(w));

            if (bound < S2[w] &&
                (!G.isDirected() || component[w] == component[source])) {
                S2[w] = bound;
            }
        }
    }
}

void DynTopHarmonicCloseness::init() {
    n = G.upperNodeIdBound();

    assert(n >= k);

    topk.clear();
    topk.resize(k);
    topkScores.clear();
    topkScores.resize(k);

    nMinCloseness = 0;
    minCloseness = std::numeric_limits<double>::max();
    trail = 0;
}

void DynTopHarmonicCloseness::run() {
    init();
    std::vector<bool> toAnalyze(n, true);
    // We compute the number of reachable nodes (or an upper bound)
    // only if we use the algorithm for complex networks.
    if (!useBFSbound) {
        computeReachableNodes();
    }

    // Main priority queue with all nodes in order of decreasing degree
    Aux::PrioQueue<edgeweight, node> Q1(n);

    G.forNodes([&](node v) { Q1.insert(-(n + G.degree(v)), v); });

    Aux::PrioQueue<edgeweight, node> top(n);

    // protects accesses to all shared variables
    omp_lock_t lock;
    omp_init_lock(&lock);

    edgeweight kth = 0;
#pragma omp parallel
    {
        std::vector<uint8_t> visited(n, false);
        std::vector<count> distances(n);
        std::vector<node> pred(n, 0);

        std::vector<edgeweight> S(n, std::numeric_limits<edgeweight>::max());

        while (!Q1.empty()) {

            omp_set_lock(&lock);
            if (Q1.empty()) { // The size of Q1 might have changed.
                omp_unset_lock(&lock);
                break;
            }
            std::pair<edgeweight, node> extracted = Q1.extractMin();

            node v = extracted.second;
            toAnalyze[v] = false;

            omp_unset_lock(&lock);

            // for networks with large diameters: break if the score of the
            // current node is smaller than the k-th highest score
            if (useBFSbound && allScores[v] < kth && isValid[v]) {
                break;
            }

            if (G.degreeOut(v) == 0) {
                omp_set_lock(&lock);
                allScores[v] = 0;
                isExact[v] = 0;
                omp_unset_lock(&lock);
            } else if (useBFSbound) {
                count visEdges;

                // Perform a complete BFS from v and obtain upper bounds
                // for all nodes in the graph
                BFSbound(v, S, &visEdges);

                omp_set_lock(&lock);
                allScores[v] = S[v];
                isExact[v] = true;
                isValid[v] = true;
                omp_unset_lock(&lock);

                // Update the scores of all nodes with the bounds obtained
                // by the complete BFS
                G.forNodes([&](node u) {
                    if (allScores[u] > S[u] &&
                        toAnalyze[u]) { // This part must be syncrhonized.
                        omp_set_lock(&lock);
                        if (allScores[u] > S[u] &&
                            toAnalyze[u]) { // Have to check again, because the variables
                            // might have changed
                            allScores[u] = S[u];
                            isValid[u] = true;
                            Q1.remove(-u);
                            Q1.insert(-allScores[u], u);
                        }
                        omp_unset_lock(&lock);
                    }
                });
            } else {
                count edgeCount = 0;

                // perform a pruned BFS from v in complex networks
                std::pair<edgeweight, bool> c =
                    BFScut(v, kth, n, r[v], visited, distances, pred, edgeCount);

                omp_set_lock(&lock);
                allScores[v] = c.first;
                isExact[v] = c.second;
                isValid[v] = true;
                omp_unset_lock(&lock);
            }

            omp_set_lock(&lock);
            // Insert v into the list with the k most central nodes if
            // its score is larger than the k-th largest value
            if (isExact[v] && allScores[v] >= kth) {
                top.insert(allScores[v], v);
                if (top.size() > k) {
                    ++trail;
                    if (allScores[v] > kth) {
                        if (nMinCloseness == trail) {
                            while (top.size() > k) {
                                top.extractMin();
                            }
                            trail = 0;
                            nMinCloseness = 1;
                            if (k > 1) {
                                double pqTop = top.peekMin(0).first;
                                double pqNext = pqTop;
                                top.forElementsWhile([&]() { return pqTop == pqNext; },
                                                     [&](double curKey, node) {
                                                         pqNext = curKey;
                                                         ++nMinCloseness;
                                                     });
                            }
                        }
                    } else {
                        ++nMinCloseness;
                    }
                } else if (allScores[v] < minCloseness) {
                    minCloseness = allScores[v];
                    nMinCloseness = 1;
                } else if (allScores[v] == minCloseness) {
                    ++nMinCloseness;
                }
            }

            // Update the k-th largest value for this thread
            if (top.size() >= k) {
                std::pair<edgeweight, node> elem = top.extractMin();
                kth = elem.first;
                top.insert(elem.first, elem.second);
                if (nMinCloseness == 1) {
                    minCloseness = kth;
                }
            }
            omp_unset_lock(&lock);
        }
    }

    if (trail > 0) {
        topk.resize(k + trail);
        topkScores.resize(k + trail);
    }

    // Store the nodes and their closeness centralities
    for (int64_t j = top.size() - 1; j >= 0; --j) {
        std::pair<edgeweight, node> elem = top.extractMin();
        topk[j] = elem.second;
        topkScores[j] = elem.first;
    }

    for (count j = 0; j < k; ++j) {
        top.insert(topkScores[j], topk[j]);
    }

    for (count i = 0; i < topk.size() - 1; ++i) {
        count toSort = 1;
        while ((i + toSort) < topk.size() &&
               topkScores[i] == topkScores[i + toSort]) {
            ++toSort;
        }

        if (toSort > 1) {
            auto begin = topk.begin() + i;
            std::sort(begin, begin + toSort);
            i += toSort - 1;
        }
    }

    hasRun = true;
}

void DynTopHarmonicCloseness::update(GraphEvent event) {

    if (event.type == GraphEvent::EDGE_ADDITION) {
        addEdge(event);
    }

    if (event.type == GraphEvent::EDGE_REMOVAL) {
        removeEdge(event);
    }
}

void DynTopHarmonicCloseness::addEdge(const GraphEvent &event) {

    n = G.upperNodeIdBound();

    node eventU = event.u;
    node eventV = event.v;

    // Compute the affected nodes
    AffectedNodes affectedNodes(G, event);
    affectedNodes.run();

    std::vector<node> uniqueAffectedNodes = affectedNodes.getNodes();
    std::vector<edgeweight> distancesFromInsertion = affectedNodes.getDistances();
    std::vector<edgeweight> improvementUpperBounds =
        affectedNodes.getImprovements();

    Aux::PrioQueue<edgeweight, node> Q1(n);

    for (node w : uniqueAffectedNodes) {

        // Add the improvement bounds for large-diameter networks
        if (useBFSbound) {
            isValid[w] = true;
            isExact[w] = false;
            allScores[w] += improvementUpperBounds[w];
        } else {
            isValid[w] = false;
        }

        Q1.insert(-allScores[w], w);
    }

    allScores[eventU] = affectedNodes.closenessU;
    isExact[eventU] = true;
    isValid[eventU] = true;
    exactCutOff[eventU] = false;
    cutOff[eventU] = std::numeric_limits<edgeweight>::max();

    if (!G.isDirected()) {
        isValid[eventV] = true;
        isExact[eventV] = true;
        allScores[eventV] = affectedNodes.closenessV;
        exactCutOff[eventV] = false;
        cutOff[eventV] = std::numeric_limits<edgeweight>::max();
    }

    Aux::PrioQueue<edgeweight, node> top(n);

    count validNodes = 0;
    trail = 0;
    nMinCloseness = 1;
    minCloseness = std::numeric_limits<double>::max();
    // insert all top k nodes with valid closeness into the queue
    for (node topNode : topk) {
        if (isValid[topNode] && isExact[topNode]) {
            top.insert(allScores[topNode], topNode);
            if (validNodes == 0) {
                minCloseness = allScores[topNode];
            } else if (minCloseness == allScores[topNode]) {
                ++nMinCloseness;
            } else {
                minCloseness = std::min(minCloseness, allScores[topNode]);
                nMinCloseness = 1;
            }
            validNodes++;
            if (validNodes > k) {
                ++trail;
            }
        }
    }

    // Recompute the number of reachable nodes for complex networks
    if (!useBFSbound) {
        // Store the old numbers
        std::copy(r.begin(), r.end(), rOld.begin());
        updateReachableNodesAfterInsertion(eventU, eventV);
    }

    // protects accesses to shared data structures
    omp_lock_t lock;
    omp_init_lock(&lock);

    // protects the variables collecting statistics
    omp_lock_t statsLock;
    omp_init_lock(&statsLock);
    // the new minimum upper bound will be larger or equal to the current
    // kth-highest closeness
    edgeweight kth = topkScores.back();
#pragma omp parallel
    {
        std::vector<uint8_t> visited(n, false);
        std::vector<count> distances(n);
        count visitedEdges = 0;
        std::vector<node> pred(n, 0);

        std::vector<edgeweight> S(n, std::numeric_limits<edgeweight>::max());

        while (!Q1.empty()) {

            omp_set_lock(&lock);
            if (Q1.empty()) {
                omp_unset_lock(&lock);
                break;
            }
            std::pair<edgeweight, node> extracted = Q1.extractMin();
            node v = extracted.second;

            omp_unset_lock(&lock);

            if (useBFSbound && allScores[v] < kth && isExact[v]) {
                break;
            }

            edgeweight distanceFromInsertedEdge = distancesFromInsertion[v];
            edgeweight boundaryUpdateScore =
                allScores[v] - 1.0 / (cutOff[v] + 2.0) + 1.0 / (cutOff[v] + 1.0);

            if ((distanceFromInsertedEdge > cutOff[v] && !isExact[v]) &&
                (!useBFSbound && r[v] <= rOld[v])) {
                // our estimate is still valid if the distance of the inserted edge is
                // larger than the previous cut-off
                isValid[v] = true;
                omp_set_lock(&statsLock);
                omp_unset_lock(&statsLock);
            } else if (distanceFromInsertedEdge == cutOff[v] && !isExact[v] &&
                       boundaryUpdateScore < kth &&
                       (!useBFSbound && r[v] <= rOld[v])) {
                // both nodes are affected but not on the same level => cheap update for
                // the upper bound
                allScores[v] = boundaryUpdateScore;
                isValid[v] = true;
                omp_set_lock(&statsLock);
                omp_unset_lock(&statsLock);
            } else if (((!useBFSbound &&
                         allScores[v] + improvementUpperBounds[v] < kth) ||
                        (useBFSbound && allScores[v] < kth))) {
                // Adding the improvement bound does not yield a higher upper bound than
                // the k-th largest value => add it
                if (!useBFSbound) {
                    allScores[v] = allScores[v] + improvementUpperBounds[v];
                    isExact[v] = false;
                    isValid[v] = true;
                }

                omp_set_lock(&statsLock);
                omp_unset_lock(&statsLock);
            } else if (!(isValid[v] && isExact[v])) {
                if (useBFSbound) {

                    if (isValid[v] && allScores[v] < kth) {
                        continue;
                    }

                    count visEdges;
                    BFSbound(v, S, &visEdges);

                    omp_set_lock(&lock);
                    allScores[v] = S[v];
                    isExact[v] = true;
                    isValid[v] = true;
                    omp_unset_lock(&lock);

                    omp_set_lock(&statsLock);
                    omp_unset_lock(&statsLock);

                    G.forNodes([&](node u) {
                        if ((allScores[u] > S[u] &&
                             !isExact[u])) { // This part must be synchronized.
                            omp_set_lock(&lock);
                            if ((allScores[u] > S[u] &&
                                 !isExact[u])) { // Have to check again, because the variables
                                // might have changed

                                allScores[u] = S[u];
                                isValid[u] = true;
                                isExact[u] = false;
                                Q1.remove(-u);
                                Q1.insert(allScores[u], u);
                            }
                            omp_unset_lock(&lock);
                        }
                    });
                } else {

                    std::pair<edgeweight, bool> c =
                        BFScut(v, kth, n, r[v], visited, distances, pred, visitedEdges);

                    allScores[v] = c.first;
                    isExact[v] = c.second;
                    isValid[v] = true;

                    omp_set_lock(&statsLock);
                    omp_unset_lock(&statsLock);
                }
            }

            omp_set_lock(&lock);
            if (isExact[v] && allScores[v] >= kth) {
                count prevSize = top.size();
                top.insert(allScores[v], v);
                if (top.size() > std::max(k, prevSize)) {
                    ++trail;
                    if (allScores[v] > kth) {
                        if (nMinCloseness == trail) {
                            while (top.size() > k) {
                                top.extractMin();
                            }
                            trail = 0;
                            nMinCloseness = 1;
                            if (k > 1) {
                                double pqTop = top.peekMin(0).first;
                                double pqNext = pqTop;
                                top.forElementsWhile([&]() { return pqTop == pqNext; },
                                                     [&](double curKey, node) {
                                                         pqNext = curKey;
                                                         ++nMinCloseness;
                                                     });
                            }
                        }
                    } else {
                        ++nMinCloseness;
                    }
                } else if (allScores[v] < minCloseness) {
                    minCloseness = allScores[v];
                    nMinCloseness = 1;
                } else if (allScores[v] == minCloseness) {
                    ++nMinCloseness;
                }
            }

            // Update the k-th largest value for this thread
            if (top.size() >= k) {
                std::pair<edgeweight, node> elem = top.extractMin();
                kth = elem.first;
                top.insert(elem.first, elem.second);
                if (nMinCloseness == 1) {
                    minCloseness = kth;
                }
            }
            omp_unset_lock(&lock);
        }
    }

    if (trail > 0) {
        topk.resize(k + trail);
        topkScores.resize(k + trail);
    }

    // Store the nodes and their closeness centralities
    for (int64_t j = top.size() - 1; j >= 0; --j) {
        std::pair<edgeweight, node> elem = top.extractMin();
        topk[j] = elem.second;
        topkScores[j] = elem.first;
    }

    for (count i = 0; i < topk.size() - 1; ++i) {
        count toSort = 1;
        while ((i + toSort) < topk.size() &&
               topkScores[i] == topkScores[i + toSort]) {
            ++toSort;
        }
        if (toSort > 1) {
            auto begin = topk.begin() + i;
            std::sort(begin, begin + toSort);
            i += toSort - 1;
        }
    }
}

void DynTopHarmonicCloseness::removeEdge(const GraphEvent &event) {

    n = G.upperNodeIdBound();

    AffectedNodes affectedNodes(G, event);
    affectedNodes.run();

    std::vector<node> uniqueAffectedNodes = affectedNodes.getNodes();

    for (node w : uniqueAffectedNodes) {
        isExact[w] = false;
    }

    Aux::PrioQueue<edgeweight, node> top(n);

    // we can abort the update early if none of the top k nodes is actually
    // affected
    bool hasAffectedTopNode = false;
    for (node v : topk) {
        if (!isExact[v]) {
            hasAffectedTopNode = true;
        }
    }

    if (!hasAffectedTopNode) {
        return;
    }

    edgeweight kth = 0;
    minCloseness = std::numeric_limits<double>::max();
    nMinCloseness = 1;
    trail = 0;

    // This yields to seg. fault when on edge removals on street networks.
    if (!useBFSbound) {
        // Since we use BFScut() here, we need the number of reachable nodes in any
        // case And now? EA
        updateReachableNodesAfterDeletion(event);
    }
    Aux::PrioQueue<edgeweight, node> Q(n);

    // add all nodes to the queue
    G.forNodes([&](node v) { Q.insert(-allScores[v], v); });

    omp_lock_t lock;
    omp_init_lock(&lock);

#pragma omp parallel
    {
        std::vector<uint8_t> visited(n, false);
        std::vector<count> distances(n);
        count visitedEdges = 0;
        std::vector<node> pred(n, 0);

        std::vector<edgeweight> S(n, std::numeric_limits<edgeweight>::max());
        while (!Q.empty()) {

            omp_set_lock(&lock);
            if (Q.empty()) {
                omp_unset_lock(&lock);
                break;
            }
            std::pair<edgeweight, node> elem = Q.extractMin();
            node v = elem.second;

            omp_unset_lock(&lock);
            if (allScores[v] < kth) {
                break;
            }

            if (!isExact[v] || !isValid[v]) {
                if (!useBFSbound) {
                    std::pair<edgeweight, bool> c =
                        BFScut(v, kth, n, r[v], visited, distances, pred, visitedEdges);

                    allScores[v] = c.first;
                    isExact[v] = c.second;
                    isValid[v] = true;
                } else {
                    BFSbound(v, S, &visitedEdges);

                    allScores[v] = S[v];
                    isExact[v] = true;
                    isValid[v] = true;
                }
            }
            omp_set_lock(&lock);
            if (isExact[v] && allScores[v] >= kth) {
                // Insert node into top queue if closeness is larger than kth
                count prevSize = top.size();
                top.insert(allScores[v], v);
                if (top.size() > std::max(k, prevSize)) {
                    ++trail;
                    if (allScores[v] > kth) {
                        if (nMinCloseness == trail) {
                            while (top.size() > k) {
                                top.extractMin();
                            }
                            trail = 0;
                            nMinCloseness = 1;
                            if (k > 1) {
                                double pqTop = top.peekMin(0).first;
                                double pqNext = pqTop;
                                top.forElementsWhile([&]() { return pqTop == pqNext; },
                                                     [&](double curKey, node) {
                                                         pqNext = curKey;
                                                         ++nMinCloseness;
                                                     });
                            }
                        }
                    } else {
                        ++nMinCloseness;
                    }
                } else if (allScores[v] < minCloseness) {
                    minCloseness = allScores[v];
                    nMinCloseness = 1;
                } else if (allScores[v] == minCloseness) {
                    ++nMinCloseness;
                }
            }

            // Update kth
            if (top.size() >= k) {
                std::pair<edgeweight, node> elem = top.extractMin();
                kth = elem.first;
                top.insert(elem.first, elem.second);
                if (nMinCloseness == 1) {
                    minCloseness = kth;
                }
            }
            omp_unset_lock(&lock);
        }
    }

    // TODO This could go to another method
    if (trail > 0) {
        topk.resize(k + trail);
        topkScores.resize(k + trail);
    }

    for (int64_t j = top.size() - 1; j >= 0; j--) {
        std::pair<edgeweight, node> elem = top.extractMin();
        topk[j] = elem.second;
        topkScores[j] = elem.first;
    }

    for (count j = 0; j < k; j++) {
        top.insert(topkScores[j], topk[j]);
    }

    for (count i = 0; i < topk.size() - 1; ++i) {
        count toSort = 1;
        while ((i + toSort) < topk.size() &&
               topkScores[i] == topkScores[i + toSort]) {
            ++toSort;
        }
        if (toSort > 1) {
            auto begin = topk.begin() + i;
            std::sort(begin, begin + toSort);
            i += toSort - 1;
        }
    }
}

void DynTopHarmonicCloseness::updateBatch(
    const std::vector<GraphEvent> &batch) {
    for (auto event : batch) {
        update(event);
    }
}

void DynTopHarmonicCloseness::reset() {
    std::fill(isValid.begin(), isValid.end(), false);
    std::fill(isExact.begin(), isExact.end(), false);
    std::fill(cutOff.begin(), cutOff.end(),
              std::numeric_limits<edgeweight>::max());
}

std::vector<std::pair<node, edgeweight>>
DynTopHarmonicCloseness::ranking(bool includeTrail) {
    count nTop = includeTrail ? k + trail : k;
    std::vector<std::pair<node, edgeweight>> ranking(nTop);

    for (count i = 0; i < nTop; ++i) {
        ranking[i] = std::make_pair(topk[i], topkScores[i]);
    }

    return ranking;
}

void DynTopHarmonicCloseness::computeReachableNodes() {
    if (G.isDirected()) {
        computeReachableNodesDirected();
    } else {
        computeReachableNodesUndirected();
    }
}

void DynTopHarmonicCloseness::computeReachableNodesUndirected() {
    if (!hasComps) {
        comps = new DynConnectedComponents(G);
        hasComps = true;
    }

    r = std::vector<count>(G.upperNodeIdBound());

    comps->run();
    std::map<index, count> sizes = comps->getComponentSizes();
    G.forNodes([&](node v) {
        index cv = comps->componentOfNode(v);
        component[v] = cv;
        r[v] = sizes[cv];
    });
}

void DynTopHarmonicCloseness::computeReachableNodesDirected() {

    r = std::vector<count>(G.upperNodeIdBound());
    wComps = new DynWeaklyConnectedComponents(G);
    wComps->run();
    hasWComps = true;
    std::map<index, count> sizes = wComps->getComponentSizes();
    G.forNodes([&](node w) {
        index cw = wComps->componentOfNode(w);
        component[w] = cw;
        r[w] = sizes[cw];
    });
}

// TODO: merge in a single method
void DynTopHarmonicCloseness::updateReachableNodesAfterInsertion(node u,
                                                                 node v) {

    GraphEvent e(GraphEvent::EDGE_ADDITION, u, v);
    if (!G.isDirected()) {

        comps->update(e);

        std::map<index, count> sizes = comps->getComponentSizes();
        DEBUG(sizes);
        G.forNodes([&](node v) {
            index cv = comps->componentOfNode(v);
            component[v] = cv;
            r[v] = sizes[cv];
        });
    } else {
        wComps->update(e);
        // TODO use alias with components so we do not have to replicate the code
        std::map<index, count> sizes = wComps->getComponentSizes();
        G.forNodes([&](node) {
            index cv = wComps->componentOfNode(v);
            component[v] = cv;
            r[v] = sizes[cv];
        });
    }
}

void DynTopHarmonicCloseness::updateReachableNodesAfterDeletion(
    const GraphEvent &event) {

    if (!G.isDirected()) {
        comps->update(event);

        std::map<index, count> sizes = comps->getComponentSizes();
        G.forNodes([&](node w) {
            index cv = comps->componentOfNode(w);
            component[w] = cv;
            r[w] = sizes[cv];
        });
    } else {
        wComps->update(event);
        std::map<index, count> sizes = wComps->getComponentSizes();
        G.forNodes([&](node w) {
            index cv = wComps->componentOfNode(w);
            component[w] = cv;
            r[w] = sizes[cv];
        });
    }
}

} /* namespace NetworKit */
