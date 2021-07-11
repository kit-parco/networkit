/*
 * TopHarmonicCloseness.cpp
 *
 * Created on: 25.02.2018
 *		 Author: nemes, Eugenio Angriman
 */

// networkit-format

#ifndef NDEBUG
#include <algorithm>
#endif
#include <atomic>
#include <omp.h>
#include <queue>

#include <networkit/centrality/TopHarmonicCloseness.hpp>
#include <networkit/reachability/ReachableNodes.hpp>

namespace NetworKit {

TopHarmonicCloseness::TopHarmonicCloseness(const Graph &G, count k, bool useNBbound)
    : G(&G), k(k), useNBbound(useNBbound), prioQ(Aux::GreaterInVector<double>{hCloseness}),
      topKNodesPQ(Aux::LessInVector<double>{hCloseness}) {

    if (k == 0 || k > G.numberOfNodes())
        throw std::runtime_error("Error: k must be in [1,...,n].");

    if (useNBbound && G.isWeighted())
        WARN("NBbound only works with unweighted graphs, edge weights will be ignored!");

    const count n = G.upperNodeIdBound();

    hCloseness.resize(n);
    reachableNodes.resize(n);

    if (!useNBbound && G.isWeighted()) {
        distanceGlobal.resize(omp_get_max_threads(), std::vector<edgeweight>(n));
        dijkstraHeaps.reserve(omp_get_max_threads());
        for (int i = 0; i < omp_get_max_threads(); ++i)
            dijkstraHeaps.emplace_back(Aux::LessInVector<edgeweight>{distanceGlobal[i]});
        minEdgeWeight = std::numeric_limits<edgeweight>::max();
        G.forEdges([&](node, node, edgeweight ew) { minEdgeWeight = std::min(minEdgeWeight, ew); });
    } else {
        visitedGlobal.resize(omp_get_max_threads(), std::vector<uint8_t>(n));
        tsGlobal.resize(omp_get_max_threads(), 0);
    }

    topKNodesPQ.reserve(k);
    prioQ.reserve(n);
    omp_init_lock(&lock);
}

TopHarmonicCloseness::~TopHarmonicCloseness() = default;

void TopHarmonicCloseness::computeReachableNodes() {
    if (G->isDirected()) {
        wccPtr = std::make_unique<WeaklyConnectedComponents>(*G);
        wccPtr->run();
        const auto compSizes = wccPtr->getComponentSizes();
        G->parallelForNodes(
            [&](node u) { reachableNodes[u] = compSizes.at(wccPtr->componentOfNode(u)); });
    } else {
        ReachableNodes rn(*G);
        rn.run();
        G->parallelForNodes([&reachableNodes = reachableNodes, &rn = rn](node u) {
            reachableNodes[u] = rn.numberOfReachableNodes(u);
        });
    }
}

void TopHarmonicCloseness::run() {
    computeReachableNodes();
    if (useNBbound)
        runNBbound();
    else
        runNBcut();

    topKNodes.resize(k + trail.size());
    topKScores.resize(k + trail.size());
    count i = k;
    do {
        --i;
        const node u = topKNodesPQ.extract_top();
        topKNodes[i] = u;
        topKScores[i] = hCloseness[u];
    } while (!topKNodesPQ.empty());

    for (; i < trail.size(); ++i) {
        topKNodes[k + i] = trail[i];
        topKScores[k + i] = hCloseness[trail[i]];
    }

    hasRun = true;
}

void TopHarmonicCloseness::runNBcut() {
    if (G->isWeighted())
        G->parallelForNodes([&](node u) { hCloseness[u] = initialBoundNBcutWeighted(u); });
    else
        G->parallelForNodes([&](node u) { hCloseness[u] = initialBoundNBcutUnweighted(u); });

    prioQ.build_heap(G->nodeRange().begin(), G->nodeRange().end());

    std::atomic_bool stop{false};
    std::atomic<double> kthCloseness{-1};

#pragma omp parallel
    while (!stop.load(std::memory_order_relaxed)) {
        node u = none;
        omp_set_lock(&lock);
        if (!prioQ.empty()) {
            u = prioQ.extract_top();
            if (topKNodesPQ.size() == k)
                kthCloseness.store(hCloseness[topKNodesPQ.top()], std::memory_order_relaxed);
            if (hCloseness[u] < kthCloseness) {
                stop.store(true, std::memory_order_relaxed);
                u = none;
            }
        } else
            stop.store(true, std::memory_order_relaxed);
        omp_unset_lock(&lock);

        if (u == none)
            break;

        if (G->isWeighted()) {
            if (!bfscutWeighted(u, kthCloseness.load(std::memory_order_relaxed)))
                continue;
        } else {
            if (!bfscutUnweighted(u, kthCloseness.load(std::memory_order_relaxed)))
                continue;
        }

        omp_set_lock(&lock);
        updateTopkPQ(u);
        omp_unset_lock(&lock);
    }
}

void TopHarmonicCloseness::runNBbound() {
    numberOfNodesAtLevelGlobal.resize(omp_get_max_threads(),
                                      std::vector<count>(G->numberOfNodes(), 0));
    nodesAtLevelGlobal.resize(omp_get_max_threads(),
                              std::vector<std::vector<count>>(G->numberOfNodes()));
    nodesAtCurrentLevelGlobal.resize(omp_get_max_threads());

    for (int i = 0; i < omp_get_max_threads(); ++i)
        nodesAtCurrentLevelGlobal[i].reserve(G->numberOfNodes() - 1);

    levelImprovement.resize(G->numberOfNodes());

    if (G->isDirected())
        computeReachableNodesBounds();
    computeNeighborhoodBasedBound();

    prioQ.build_heap(G->nodeRange().begin(), G->nodeRange().end());

    std::atomic_bool stop{false};

#pragma omp parallel
    while (!stop.load(std::memory_order_relaxed)) {
        node u = none;
        omp_set_lock(&lock);
        if (!prioQ.empty()) {
            u = prioQ.extract_top();
            if (topKNodes.size() == k && hCloseness[u] <= hCloseness[topKNodesPQ.top()]) {
                stop.store(true, std::memory_order_relaxed);
                u = none;
            }
        } else
            stop.store(true, std::memory_order_relaxed);
        omp_unset_lock(&lock);

        if (u == none)
            break;

        bfsbound(u);

        omp_set_lock(&lock);
        updateTopkPQ(u);
        omp_unset_lock(&lock);
    }
}

void TopHarmonicCloseness::updateTopkPQ(node u) {
    topKNodesPQ.push(u);
    if (topKNodesPQ.size() <= k)
        return;

    assert(topKNodesPQ.size() == k + 1);
    const node bottomNode = topKNodesPQ.extract_top();
    const double kthCloseness = hCloseness[topKNodesPQ.top()];
    if (kthCloseness == hCloseness[bottomNode]) {
        if (trail.empty() || hCloseness[bottomNode] == hCloseness[trail[0]]) {
            trail.push_back(bottomNode);
        } else {
            assert(hCloseness[bottomNode] > hCloseness[trail[0]]);
            trail.clear();
            trail.push_back(bottomNode);
        }
    } else
        trail.clear();
}

bool TopHarmonicCloseness::bfscutUnweighted(node source, double kthCloseness) {
    const count reachableFromSource = reachableNodes[source];
    const count undirected = !G->isDirected();
    updateTimestamp();
    auto &visited = visitedGlobal[omp_get_thread_num()];
    const auto ts = tsGlobal[omp_get_thread_num()];
    visited[source] = ts;
    count visitedNodes = 1, level = 1;

    std::queue<node> q1, q2;
    q1.push(source);

    double h = 0, htilde = 0;

    do {
        count nodesAtNextLevelUB = 0;
        do {
            const node u = q1.front();
            q1.pop();

            G->forNeighborsOf(u, [&](node v) {
                if (visited[v] == ts)
                    return;
                visited[v] = ts;
                q2.push(v);
                ++visitedNodes;
                h += 1. / static_cast<double>(level);
                nodesAtNextLevelUB += G->degree(v) - undirected;
            });
        } while (!q1.empty());

        assert(reachableFromSource >= visitedNodes);
        nodesAtNextLevelUB = std::min(nodesAtNextLevelUB, reachableFromSource - visitedNodes);
        htilde = h;
        htilde += static_cast<double>(nodesAtNextLevelUB) / static_cast<double>(level + 1);
        htilde += static_cast<double>(reachableFromSource - visitedNodes - nodesAtNextLevelUB)
                  / static_cast<double>(level + 2);

#ifndef NDEBUG
        if (q2.empty()) { // Finished the BFS
            if (G->isDirected())
                assert(reachableFromSource >= visitedNodes);
            else
                assert(reachableFromSource == visitedNodes);
        }
#endif

        // Prune BFS
        if (htilde < kthCloseness) {
            hCloseness[source] = htilde;
            return false;
        }

        ++level;
        std::swap(q1, q2);
    } while (!q1.empty());

    hCloseness[source] = h;
    return true;
}

bool TopHarmonicCloseness::bfscutWeighted(node source, double kthCloseness) {
    auto &distance = distanceGlobal[omp_get_thread_num()];
    std::fill(distance.begin(), distance.end(), std::numeric_limits<edgeweight>::max());
    distance[source] = 0;
    auto &pq = dijkstraHeaps[omp_get_thread_num()];
    pq.clear();
    const count reachableFromSource = reachableNodes[source];

    // Visit the source immediately so we can avoid 'if (u != source) ... ' in the while loop below
    G->forNeighborsOf(source, [&](node v, edgeweight ew) {
        distance[v] = ew;
        pq.update(v);
    });

    double h = 0, htilde = 0;
    count visitedNodes = 1;

    while (!pq.empty()) {
        const node u = pq.extract_top();
        assert(u != source);
        assert(distance[u] > 0);

        ++visitedNodes;
        assert(reachableFromSource >= visitedNodes);
        const double distU = distance[u];
        h += 1. / distU;
        // Simple upper bound for harmonic closeness: we assume that the distance from 'source' to
        // all the remaining unvisited nodes is d(source, u).
        htilde = h + static_cast<double>(reachableFromSource - visitedNodes) / distU;

        // Prune SSSP
        if (htilde < kthCloseness) {
            hCloseness[source] = htilde;
            return false;
        }

        G->forNeighborsOf(u, [&](node v, edgeweight ew) {
            const double newDistV = distU + ew;
            if (newDistV < distance[v]) {
                distance[v] = newDistV;
                pq.update(v);
            }
        });
    }

#ifndef NDEBUG
    if (G->isDirected())
        assert(reachableFromSource >= visitedNodes);
    else
        assert(reachableFromSource == visitedNodes);
#endif

    hCloseness[source] = h;
    return true;
}

void TopHarmonicCloseness::bfsbound(node source) {
    updateTimestamp();
    auto &visited = visitedGlobal[omp_get_thread_num()];
    const auto ts = tsGlobal[omp_get_thread_num()];
    visited[source] = ts;
    auto &nodesAtLevel = nodesAtLevelGlobal[omp_get_thread_num()];
    nodesAtLevel.clear();
    auto &nodesAtCurrentLevel = nodesAtCurrentLevelGlobal[omp_get_thread_num()];
    nodesAtCurrentLevel.clear();
    auto &numberOfNodesAtLevel = numberOfNodesAtLevelGlobal[omp_get_thread_num()];

    for (count i = 1; i < numberOfNodesAtLevel.size() && numberOfNodesAtLevel[i] > 0; ++i)
        numberOfNodesAtLevel[i] = 0;

    assert(std::count_if(numberOfNodesAtLevel.begin(), numberOfNodesAtLevel.end(),
                         [](count i) { return i > 0; })
           == 0);

    std::queue<node> q1, q2;
    q1.push(source);
    count level = 1;
    hCloseness[source] = 0;

    do {
        do {
            const node u = q1.front();
            q1.pop();
            nodesAtCurrentLevel.push_back(u);

            G->forNeighborsOf(u, [&](node v) {
                if (visited[v] == ts)
                    return;
                visited[v] = ts;
                q2.push(v);
                ++numberOfNodesAtLevel[level];
                hCloseness[source] += 1. / static_cast<double>(level);
            });
        } while (!q1.empty());

        nodesAtLevel.emplace_back(std::move(nodesAtCurrentLevel));
        nodesAtCurrentLevel.clear();
        ++level;
        std::swap(q1, q2);
    } while (!q1.empty());

    const double nearNodes = 1. + numberOfNodesAtLevel[1] + numberOfNodesAtLevel[2];
    double farNodes = 0;
    for (count i = 3; i < level; ++i)
        farNodes += static_cast<double>(numberOfNodesAtLevel[i]) / static_cast<double>(i - 1);

    double levelBound = nearNodes / 2. + farNodes;

    const auto updateBound = [&](const node y) -> void {
        // y has already been processed
        if (!prioQ.contains(y))
            return;
        const double bound = levelBound + (static_cast<double>(G->degree(y)) - 1.) / 2.;
        if (bound < hCloseness[y]) {
            if (G->isDirected()) {
                if (wccPtr->componentOfNode(source) == wccPtr->componentOfNode(y)) {
                    hCloseness[y] = bound;
                    prioQ.update(y);
                }
            } else {
                hCloseness[y] = bound;
                prioQ.update(y);
            }
        }
    };

    omp_set_lock(&lock);
    // Update bound for vertices at level 1
    G->forNeighborsOf(source, updateBound);
    omp_unset_lock(&lock);

    // Update bound for vertices at level >= 2
    for (int64_t i = 2; i < static_cast<int64_t>(level); ++i) {
        for (int64_t j = 0; j < static_cast<int64_t>(level); ++j) {
            double denominator = 2;
            if (G->isDirected())
                denominator = j - i > 2 ? static_cast<double>(j - i) : denominator;
            else
                denominator = std::max(denominator, static_cast<double>(std::abs(j - i)));
            levelBound += static_cast<double>(numberOfNodesAtLevel[j]) / denominator;
        }

        if (static_cast<size_t>(i) <= nodesAtLevel.size()) {
            omp_set_lock(&lock);
            for (node y : nodesAtLevel[i - 1])
                updateBound(y);
            omp_unset_lock(&lock);
        }
    }
}

double TopHarmonicCloseness::initialBoundNBcutUnweighted(node u) const noexcept {
    const count degU = G->degree(u);
    if (degU == 0)
        return 0;
    return static_cast<double>(degU) + static_cast<double>(reachableNodes[u] - degU) / 2.;
}

double TopHarmonicCloseness::initialBoundNBcutWeighted(node u) const noexcept {
    if (G->degree(u) == 0)
        return 0;
    edgeweight minOutEdgeWeight = std::numeric_limits<edgeweight>::max();
    G->forNeighborsOf(
        u, [&](node, edgeweight ew) { minOutEdgeWeight = std::min(minOutEdgeWeight, ew); });
    return 1. / minOutEdgeWeight
           + static_cast<edgeweight>(reachableNodes[u] - 1) / (minOutEdgeWeight + minEdgeWeight);
}

template <class Type>
std::vector<Type> TopHarmonicCloseness::resizedCopy(const std::vector<Type> &vec) const noexcept {
    std::vector<Type> result = vec;
    result.resize(k);
    return result;
}

std::vector<node> TopHarmonicCloseness::topkNodesList(bool includeTrail) {
    if (!includeTrail && !trail.empty())
        return resizedCopy(topKNodes);
    return topKNodes;
}

std::vector<edgeweight> TopHarmonicCloseness::topkScoresList(bool includeTrail) {
    if (!includeTrail && !trail.empty())
        return resizedCopy(topKScores);
    return topKScores;
}

void TopHarmonicCloseness::computeReachableNodesBounds() {
    ReachableNodes rn(*G, false);
    rn.run();

    reachL.resize(G->upperNodeIdBound());
    reachU.resize(G->upperNodeIdBound());

    G->parallelForNodes([&rn = rn, &reachL = reachL, &reachU = reachU](node u) {
        reachL[u] = rn.numberOfReachableNodesLB(u);
        reachU[u] = rn.numberOfReachableNodesUB(u);
    });
}

void TopHarmonicCloseness::computeNeighborhoodBasedBound() {
    const count n = G->upperNodeIdBound();
    std::vector<count> nodesAtK(n);
    std::vector<count> nodesAtKminusOne(n);
    std::vector<count> nodesAtKminusTwo(n);
    std::vector<count> visited(n);
    std::vector<bool> finished(n);
    std::vector<double> sumHCloseness(n);

    count nFinished = 0;
#pragma omp parallel for reduction(+ : nFinished)
    for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
        const node u = static_cast<node>(i);
        if (!G->hasNode(u))
            continue;
        const count degOutU = G->degreeOut(u);
        if (degOutU == 0) {
            ++nFinished;
            finished[u] = true;
        }

        nodesAtKminusOne[u] = degOutU;
        sumHCloseness[u] = degOutU;
        visited[u] = degOutU + 1;
        hCloseness[u] = std::numeric_limits<double>::max();
    }

    count level = 2;
    while (nFinished < G->numberOfNodes()) {
        G->forNodes([&](const node u) {
            if (finished[u])
                return;

            nodesAtK[u] = 0;
            G->forNeighborsOf(u, [&](node v) { nodesAtK[u] += nodesAtKminusOne[v]; });
            if (!G->isDirected()) {
                if (level == 2) {
                    assert(nodesAtK[u] >= G->degree(u));
                    nodesAtK[u] -= G->degree(u);
                } else {
                    assert(nodesAtK[u] >= nodesAtKminusTwo[u] * (G->degreeOut(u) - 1));
                    nodesAtK[u] -= nodesAtKminusTwo[u] * (G->degreeOut(u) - 1);
                }
            }

            const count nOld = visited[u];
            visited[u] += nodesAtK[u];
            sumHCloseness[u] += static_cast<double>(nodesAtK[u]) / static_cast<double>(level);

            if (G->isDirected()) {
                if (visited[u] >= reachL[u]) {
                    if (nOld < reachL[u]) {
                        // We have to consider the case in which the number of reachable
                        // vertices is reachL.
                        hCloseness[u] = sumHCloseness[u]
                                        - static_cast<double>(visited[u] - reachL[u])
                                              / static_cast<double>(level);
                    }

                    if (nodesAtK[u] == 0)
                        reachU[u] = visited[u];

                    if (visited[u] >= reachU[u]) {
                        // We have to consider the case in which the number of reachable
                        // vertices is reachU.
                        hCloseness[u] = std::min(hCloseness[u],
                                                 sumHCloseness[u]
                                                     - static_cast<double>(visited[u] - reachU[u])
                                                           / static_cast<double>(level));
                        finished[u] = true;
                        nFinished++;

                        assert(visited[u] >= reachL[u] || nodesAtK[u] != 0);
                    } else { // reachL < N < reachU
                        // We have to consider the case in which the number of reachable is
                        // N[u].
                        hCloseness[u] = std::min(hCloseness[u], sumHCloseness[u]);
                    }
                }
            } else if (visited[u] >= reachableNodes[u]) {
                hCloseness[u] = std::min(hCloseness[u],
                                         sumHCloseness[u]
                                             - static_cast<double>(visited[u] - reachableNodes[u])
                                                   / static_cast<double>(level));
                finished[u] = true;
                ++nFinished;
            }
        });

        G->parallelForNodes([&](node u) {
            nodesAtKminusTwo[u] = nodesAtKminusOne[u];
            nodesAtKminusOne[u] = nodesAtK[u];
        });
        ++level;
    }
}

} // namespace NetworKit
