/*
 * TopCloseness.cpp
 *
 *  Created on: 03.10.2014
 *      Author: ebergamini, michele borassi
 */

// networkit-format

#include <omp.h>
#include <queue>
#include <stack>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/VectorComparator.hpp>
#include <networkit/centrality/TopCloseness.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/reachability/ReachableNodes.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>
#include <tlx/container/d_ary_heap.hpp>

namespace NetworKit {

TopCloseness::TopCloseness(const Graph &G, count k, bool first_heu, bool sec_heu)
    : G(G), k(k), first_heu(first_heu), sec_heu(sec_heu) {}

void TopCloseness::init() {
    n = G.upperNodeIdBound();
    trail = 0;
    visEdges = 0;
    n_op = 0;
    maxFarness = 0.0;
    nMaxFarness = 0;
    topk.clear();
    topk.resize(k);
    topkScores.clear();
    topkScores.resize(k);
    DEBUG("Number of nodes: ", n);
    DEBUG("k = ", k);
    farness.clear();
    farness.resize(n, 0);
    if (sec_heu) {
        nodesPerLevs.resize(omp_get_max_threads(), std::vector<count>(n));
        sumLevels.resize(omp_get_max_threads(), std::vector<count>(n));
    }
    computeReachable();
    if (G.isDirected()) {
        sccsPtr = std::make_unique<StronglyConnectedComponents>(G);
        sccsPtr->run();
    }
    DEBUG("Done INIT");
}

void TopCloseness::computeReachable() {
    ReachableNodes rn(G, false);
    rn.run();

    reachLPtr = std::make_shared<std::vector<count>>(n);
    reachUPtr = std::make_shared<std::vector<count>>(n);
    auto &reachL = *reachLPtr, &reachU = *reachUPtr;

    G.parallelForNodes([&](node u) {
        reachL[u] = rn.numberOfReachableNodesLB(u);
        reachU[u] = rn.numberOfReachableNodesUB(u);
    });
}

void TopCloseness::computelBound1(std::vector<double> &S) {
    std::vector<count> neighbors(n, 0);
    std::vector<count> N(n, 0);
    std::vector<count> neighbors_new(n, 0);
    std::vector<count> neighbors_old(n, 0);
    std::vector<double> sumDist(n, 0);
    std::vector<bool> finished(n, false);

    count n_finished = 0;

    G.forNodes([&](node u) {
        S[u] = std::numeric_limits<double>::max();
        if (G.degreeOut(u) == 0) {
            finished[u] = true;
            n_finished++;
        }
        neighbors[u] = G.degreeOut(u);
        sumDist[u] = neighbors[u];
        N[u] = neighbors[u] + 1; // we also count the node itself in the number of visited nodes
    });
    count level = 2;
    DEBUG("computing first lbound");

    auto &reachL = *reachLPtr, &reachU = *reachUPtr;
    while (n_finished < n) {
        DEBUG("First bound. Finished: ", n_finished, " of ", n, ".");
        G.forNodes([&](node u) {
            if (!finished[u]) {
                n_op += G.degreeOut(u);
                neighbors_new[u] = 0;
                G.forNeighborsOf(u, [&](node v) { neighbors_new[u] += neighbors[v]; });
                if (!G.isDirected()) {
                    if (level == 2) {
                        neighbors_new[u] -= G.degreeOut(u);
                    } else {
                        if (neighbors_new[u] < (G.degreeOut(u) - 1) * neighbors_old[u]) {
                            DEBUG("BIG MISTAKE");
                            while (true) {
                            }
                        }
                        neighbors_new[u] -= (G.degreeOut(u) - 1) * neighbors_old[u];
                    }
                }

                count n_old = N[u];
                N[u] += neighbors_new[u];
                sumDist[u] += static_cast<double>(level * neighbors_new[u]);

                if (N[u] >= reachL[u]) {
                    if (n_old < reachL[u]) {
                        // We have to consider the case in which the number of reachable
                        // vertices is reachL.
                        double s1 = sumDist[u] - static_cast<double>(level * (N[u] - reachL[u]));
                        double s2 = static_cast<double>(n - 1) / static_cast<double>(reachL[u] - 1)
                                    / static_cast<double>(reachL[u] - 1);
                        S[u] = s1 * s2;
                    }
                    if (neighbors_new[u] == 0) {
                        reachU[u] = N[u];
                    }
                    if (N[u] >= reachU[u]) {
                        // We have to consider the case in which the number of reachable
                        // vertices is reachU.
                        S[u] = std::min(
                            S[u], (sumDist[u] - static_cast<double>(level * (N[u] - reachU[u])))
                                      * static_cast<double>(n - 1)
                                      / static_cast<double>(reachU[u] - 1)
                                      / static_cast<double>(reachU[u] - 1));
                        finished[u] = true;
                        n_finished++;

                        assert(N[u] >= reachL[u] || neighbors_new[u] != 0);
                    } else { // reachL < N < reachU
                        // We have to consider the case in which the number of reachable is
                        // N[u].
                        S[u] = std::min(S[u], sumDist[u] * static_cast<double>(n - 1)
                                                  / static_cast<double>(N[u] - 1)
                                                  / static_cast<double>(N[u] - 1));
                    }
                }
            }
        });
        G.forNodes([&](node u) {
            // We update neighbors.
            neighbors_old[u] = neighbors[u];
            neighbors[u] = neighbors_new[u];
        });
        level++;
    }
    DEBUG("Visited edges (first lbound): ", n_op);
}

void TopCloseness::BFSbound(node x, std::vector<double> &S2, count &visEdges,
                            const std::vector<bool> &toAnalyze) {
    count r = 0;
    std::vector<std::vector<node>> levels(n);
    // nodesPerLev[i] contains the number of nodes in level i
    std::vector<count> &nodesPerLev = nodesPerLevs[omp_get_thread_num()];
    std::fill(nodesPerLev.begin(), nodesPerLev.end(), 0);
    // sumLevs[i] contains the sum of the nodes in levels j <= i
    std::vector<count> &sumLevs = sumLevels[omp_get_thread_num()];
    std::fill(sumLevs.begin(), sumLevs.end(), 0);
    count nLevs = 0;
    double sum_dist = 0;
    Traversal::BFSfrom(G, x, [&](node u, count dist) {
        sum_dist += dist;
        r++;
        if (dist > nLevs) {
            sumLevs[nLevs] += nodesPerLev[nLevs];
            sumLevs[nLevs + 1] = sumLevs[nLevs];
            nLevs++;
            levels[nLevs].clear();
        }
        levels[nLevs].push_back(u);
        nodesPerLev[nLevs]++;
    });
    sumLevs[nLevs] += nodesPerLev[nLevs];
    if (G.isDirected()) {
        visEdges += G.numberOfEdges();
    } else {
        visEdges += 2 * G.numberOfEdges();
    }
    S2[x] = sum_dist * (n - 1.0) / (r - 1.0) / (r - 1.0);
    // we compute the bound for the first level
    count closeNodes = 0, farNodes = 0;
    for (count j = 0; j <= nLevs; j++) {
        if (std::abs(static_cast<long long>(j) - 1LL) <= 1) {
            closeNodes += nodesPerLev[j];
        } else {
            farNodes += nodesPerLev[j] * std::abs(1LL - static_cast<long long>(j));
        }
    }

    edgeweight level_bound = 2.0 * closeNodes + static_cast<double>(farNodes);
    const auto &reachU = *(reachUPtr.get());
    for (count j = 0; j < levels[1].size(); j++) {
        node w = levels[1][j];
        // we subtract 2 not to count the node itself
        double bound =
            (level_bound - 2 - G.degree(w)) * (n - 1.0) / (reachU[w] - 1.0) / (reachU[w] - 1.0);
        if (toAnalyze[w] && bound > S2[w]
            && (!G.isDirected() || sccsPtr->componentOfNode(w) == sccsPtr->componentOfNode(x))) {
            S2[w] = bound;
        }
    }

    // now we compute it for the other levels
    for (omp_index i = 2; i <= static_cast<omp_index>(nLevs); i++) {
        if (!G.isDirected() && i > 2) {
            level_bound += sumLevs[i - 3];
        }
        if (i < nLevs) {
            level_bound -= static_cast<edgeweight>(sumLevs[nLevs] - sumLevs[i + 1]);
        }
        for (count j = 0; j < levels[i].size(); j++) {
            node w = levels[i][j];
            double bound =
                (level_bound - 2 - G.degree(w)) * (n - 1.0) / (reachU[w] - 1.0) / (reachU[w] - 1.0);
            if (toAnalyze[w] && bound > S2[w]
                && (!G.isDirected()
                    || sccsPtr->componentOfNode(w) == sccsPtr->componentOfNode(x))) {
                // TODO MICHELE: as before.
                S2[w] = bound;
            }
        }
    }
}

double TopCloseness::BFScut(node v, double x, std::vector<bool> &visited,
                            std::vector<count> &distances, std::vector<node> &pred,
                            count &visEdges) {
    count d = 0, f = 0, nd = 1;
    const double rL = (*reachLPtr)[v], rU = (*reachUPtr)[v];
    std::queue<node> Q1;
    std::queue<node> to_reset;
    count sum_dist = 0;
    double ftildeL = 0, ftildeU = 0, gamma = G.degreeOut(v);
    double farnessV = 0;

    // MICHELE: variable visited is not local, otherwise the allocation would be
    // too expensive.
    visited[v] = true;
    distances[v] = 0;
    Q1.push(v);
    to_reset.push(v);

    do {
        node u = Q1.front();
        Q1.pop();

        sum_dist += distances[u];
        if (distances[u] > d) { // Need to update bounds!
            d++;
            double f1 = f + static_cast<double>(d + 2) * static_cast<double>(rL - nd) - gamma;
            double f2 = static_cast<double>(n - 1) / (rL - 1.0) / (rL - 1.0);
            ftildeL = f1 * f2;
            f1 = f + static_cast<double>(d + 2) * static_cast<double>(rU - nd) - gamma;
            f2 = static_cast<double>(n - 1) / (rU - 1.0) / (rU - 1.0);
            ftildeU = f1 * f2;
            if (std::min(ftildeL, ftildeU) >= x) {
                farnessV = std::min(ftildeL, ftildeU);
                break;
            }
            gamma = 0;
        }
        bool cont = true;
        G.forNeighborsOf(u, [&](node w) {
            if (cont) {
                ++visEdges;
                if (!visited[w]) {
                    distances[w] = distances[u] + 1;
                    Q1.push(w);
                    to_reset.push(w);
                    visited[w] = true;
                    f += distances[w];
                    if (!G.isDirected())
                        gamma += static_cast<double>(G.degree(w)
                                                     - 1); // notice: only because it's undirected
                    else
                        gamma += G.degreeOut(w);
                    nd++;
                    pred[w] = u;
                } else {
                    if (G.isDirected() || pred[u] != w) {
                        ftildeL += static_cast<double>(n - 1) / (rL - 1.0) / (rL - 1.0);
                        ftildeU += static_cast<double>(n - 1) / (rU - 1.0) / (rU - 1.0);
                        if (std::min(ftildeL, ftildeU) > x) {
                            cont = false;
                        }
                    }
                }
            }
        });
        if (std::min(ftildeL, ftildeU) > x) {
            farnessV = std::min(ftildeL, ftildeU);
            break;
        }
    } while (!Q1.empty());

    do {
        // MICHELE: need to reset variable visited.
        // Variables pred and distances
        // do not need to be updated.
        node u = to_reset.front();
        to_reset.pop();
        visited[u] = false;
    } while (!to_reset.empty());
    if (farnessV < x) {
        farnessV = static_cast<double>(sum_dist * (n - 1)) / static_cast<double>(nd - 1.0)
                   / static_cast<double>(nd - 1.0);
    }
    return farnessV;
}

void TopCloseness::run() {
    init();
    tlx::d_ary_heap<node, 2, Aux::GreaterInVector<double>> top{farness};
    top.reserve(k);
    std::vector<bool> toAnalyze(n, true);
    omp_lock_t lock;
    omp_init_lock(&lock);

    std::vector<double> S(n);
    // first lower bound on s
    if (first_heu) {
        DEBUG("Computing Neighborhood-based lower bound");
        computelBound1(S);
    }

    DEBUG("Initializing queue");
    G.forNodes([&](node u) {
        if (G.degreeOut(u) == 0) {
            farness[u] = std::numeric_limits<double>::max();
        } else if (first_heu) {
            farness[u] = S[u];
        } else {
            farness[u] = -(static_cast<double>(G.degreeOut(u)));
        }
    });
    tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<double>> Q{farness};

    std::vector<node> nodes;
    nodes.reserve(G.numberOfNodes());
    G.forNodes([&](const node u) { nodes.emplace_back(u); });
    Q.build_heap(std::move(nodes));
    DEBUG("Done filling the queue");

    std::vector<std::vector<bool>> visitedVec;
    std::vector<std::vector<count>> distVec;
    std::vector<std::vector<node>> predVec;

    if (!sec_heu) {
        visitedVec.resize(omp_get_max_threads(), std::vector<bool>(n));
        distVec.resize(omp_get_max_threads(), std::vector<count>(n));
        predVec.resize(omp_get_max_threads(), std::vector<node>(n));
    }

    double kth = std::numeric_limits<double>::max(); // like in Crescenzi
#pragma omp parallel                                 // Shared variables:
    // cc: synchronized write, read leads to a positive race condition;
    // Q: fully synchronized;
    // top: fully synchronized;
    // toAnalyze: fully synchronized;
    // visEdges: one variable for each thread, summed at the end;
    {
        count visEdges = 0;
#ifndef NETWORKIT_RELEASE_LOGGING
        count iters = 0;
#endif

        while (!Q.empty()) {
            DEBUG("To be analyzed: ", Q.size());
            omp_set_lock(&lock);
            if (Q.empty()) { // The size of Q might have changed.
                omp_unset_lock(&lock);
                break;
            }
            // Access to Q must be synchronized
            node s = Q.extract_top();
            toAnalyze[s] = false;
            omp_unset_lock(&lock);

            if (G.degreeOut(s) == 0 || farness[s] > kth) {
                break;
            }
            DEBUG("Iteration ", ++iters, " of thread ", omp_get_thread_num());

            DEBUG("    Extracted node ", s, " with priority ", farness[s], ".");
            if (G.degreeOut(s) == 0) {

                omp_set_lock(&lock);
                toAnalyze[s] = false;
                farness[s] = std::numeric_limits<double>::max();
                omp_unset_lock(&lock);

            } else if (sec_heu) {
                // MICHELE: we use BFSbound to bound the centrality of all nodes.
                DEBUG("    Running BFSbound.");
                BFSbound(s, S, visEdges, toAnalyze);
                omp_set_lock(&lock);
                farness[s] = S[s];
                omp_unset_lock(&lock);
                count imp = 0;
                for (count v = 0; v < n; v++) {
                    if (farness[v] < S[v] && toAnalyze[v]) { // This part must be syncrhonized.
                        omp_set_lock(&lock);
                        if (farness[v] < S[v] && toAnalyze[v]) { // Have to check again, because the
                                                                 // variables might have changed
                            imp++;
                            farness[v] = S[v];
                            Q.update(v);
                        }
                        omp_unset_lock(&lock);
                    }
                }
                DEBUG("    We have improved ", imp, " bounds.");
            } else {
                // MICHELE: we use BFScut to bound the centrality of s.
                DEBUG("    Running BFScut with x=", kth, " (degree:", G.degreeOut(s), ").");
                auto &visited = visitedVec[omp_get_thread_num()];
                std::fill(visited.begin(), visited.end(), false);
                auto &distances = distVec[omp_get_thread_num()];
                auto &pred = predVec[omp_get_thread_num()];
                const double farnessS = BFScut(s, kth, visited, distances, pred, visEdges);
                DEBUG("    Visited edges: ", visEdges, ".");
                omp_set_lock(&lock);
                farness[s] = farnessS;
                omp_unset_lock(&lock);
            }

            // If necessary, we update kth.
            omp_set_lock(&lock);
            if (farness[s] <= kth) {
                DEBUG("    The closeness of s is ", 1.0 / farness[s], ".");
                top.push(s);
                if (top.size() > k) {
                    ++trail;
                    if (farness[s] < kth) {
                        if (nMaxFarness == trail) {
                            // Purging trail
                            do {
                                top.extract_top();
                            } while (top.size() > k);

                            trail = 0;
                            nMaxFarness = 1;
                            if (k > 1) {
                                node last = top.extract_top();
                                maxFarness = farness[last];
                                std::stack<node> tmp;
                                tmp.push(last);

                                while (!top.empty() && farness[last] == farness[top.top()]) {
                                    tmp.push(top.extract_top());
                                    ++nMaxFarness;
                                }

                                while (!tmp.empty()) {
                                    top.push(tmp.top());
                                    tmp.pop();
                                }
                            }
                        }
                    } else { // Same farness as kth
                        ++nMaxFarness;
                    }
                } else if (farness[s] > maxFarness) {
                    maxFarness = farness[s];
                    nMaxFarness = 1;
                } else if (farness[s] == maxFarness) {
                    ++nMaxFarness;
                }
            } else {
                DEBUG("    Not in the top-k.");
            }

            // We load the new value of kth.
            if (top.size() >= k) {
                kth = farness[top.top()];
                if (nMaxFarness == 1) {
                    maxFarness = kth;
                }
            }
            omp_unset_lock(&lock);
        }
        DEBUG("Number of iterations of thread ", omp_get_thread_num(), ": ", iters, " out of ", n);
        omp_set_lock(&lock);
        this->visEdges += visEdges;
        omp_unset_lock(&lock);
    }

    if (trail) {
        topk.resize(k + trail);
        topkScores.resize(k + trail);
    }

    for (count i = top.size(); i; --i) {
        node elem = top.extract_top();
        topk[i - 1] = elem;
        topkScores[i - 1] = 1.0 / farness[elem];
    }
    for (count j = 0; j < k + trail; j++) {
        DEBUG(j + 1, "-th node with max closeness: ", topk[j], ", its closeness: ", topkScores[j]);
    }

    // Ascending order of nodes ids with same closeness
    for (count i = 0; i < topk.size() - 1; ++i) {
        count toSort = 1;
        while ((i + toSort) < topk.size() && topkScores[i] == topkScores[i + toSort]) {
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

} /* namespace NetworKit */
