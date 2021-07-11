/*
 * GroupHarmonicCloseness.cpp
 *
 * Created on: 15.12.2020
 *     Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#include <atomic>
#include <cassert>
#include <cmath>
#include <limits>
#include <omp.h>
#include <queue>
#include <stdexcept>

#include <networkit/auxiliary/VectorComparator.hpp>
#include <networkit/centrality/GroupHarmonicCloseness.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/Dijkstra.hpp>
#include <networkit/reachability/ReachableNodes.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

namespace {

struct SSSPResult {
    bool complete;
    double score;
    SSSPResult(bool complete, double score) : complete(complete), score(score) {}
};

template <class WeightType>
class GroupHarmonicClosenessImpl final
    : public GroupHarmonicCloseness::GroupHarmonicClosenessInterface {

    static constexpr WeightType infDist = std::numeric_limits<WeightType>::max();

public:
    GroupHarmonicClosenessImpl(const Graph &G, count k = 1);

    ~GroupHarmonicClosenessImpl() override;

    void run() override;

    static double scoreOfGroup(const Graph &graph, const std::vector<node> &inputGroup);

private:
    const Graph *G;
    const count k;
    WeightType minEdgeWeight = infDist;

    // Connected components partitioning
    Partition component;

    // Distance between the group to every node
    std::vector<WeightType> distFromGroup;
    // Distance between the group to every node after the most promising swap is performed
    std::vector<WeightType> curBestDist;
    // Distance vectors used by SSSPs to avoid dynamic memory allocations
    std::vector<std::vector<WeightType>> curDistGlobal;
    // Number of nodes reachable in each connected component
    std::vector<count> reachableNodesInComponent;
    // Upper bound of nodes reachable from each node
    std::vector<count> reachableNodesUB;
    // Predecessor node in a SSSP tree
    std::vector<node> pred;
    // Nodes closer to a vertex than to the group
    std::vector<node> nearerNodes;
    // As above, thread-local vector
    std::vector<std::vector<node>> curNearerNodesGlobal;
    // Marginal gain of each node
    std::vector<double> margGain;
    // Visited nodes, used in graph explorations, one vector per thread
    std::vector<std::vector<bool>> visitedGlobal;

    tlx::d_ary_addressable_int_heap<node, 2, Aux::GreaterInVector<double>> candidateNodesPQ;
    std::vector<tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<WeightType>>>
        dijkstraHeaps;

    void computeReachableNodesUndirected();
    void computeReachableNodesDirected();

    SSSPResult prunedSSSPEmptyGroup(node source, double bestScore);
    SSSPResult prunedSSSP(node source, double bestScore);

    node findTopHarmonicCloseness();
    node findNodeWithHighestMargGain();

    double harmonicClosenessUBUndirected(node u) const noexcept;
    double harmonicClosenessUBDirected(node u) const noexcept;

#ifdef NETWORKIT_SANITY_CHECKS
    void checkDistFromGroup() const;
#endif
};

template <class WeightType>
constexpr WeightType GroupHarmonicClosenessImpl<WeightType>::infDist;

template <class WeightType>
GroupHarmonicClosenessImpl<WeightType>::GroupHarmonicClosenessImpl(const Graph &G, count k)
    : G(&G), k(k), candidateNodesPQ{Aux::GreaterInVector<double>(margGain)} {

    if (k == 0 || G.numberOfNodes() <= k)
        throw std::runtime_error("Error, k must be in [1, n - 1]");

    const count n = G.upperNodeIdBound();
    const omp_index t = omp_get_max_threads();

    distFromGroup.resize(n, infDist);
    curBestDist.resize(n);
    curDistGlobal.resize(t, std::vector<WeightType>(n, infDist));
    group.reserve(k);
    nearerNodes.reserve(n - 1);
    margGain.resize(n);
    visitedGlobal.resize(t, std::vector<bool>(n, false));

    curNearerNodesGlobal.resize(t);
    for (omp_index i = 0; i < t; ++i)
        curNearerNodesGlobal[i].reserve(n - 1);

    candidateNodesPQ.reserve(n);
    if (G.isWeighted()) {
        dijkstraHeaps.reserve(t);
        for (omp_index i = 0; i < t; ++i) {
            dijkstraHeaps.emplace_back(curDistGlobal[i]);
            dijkstraHeaps.back().reserve(n);
        }
        G.forEdges([&minEdgeWeight = minEdgeWeight](node, node, WeightType ew) {
            minEdgeWeight = std::min(minEdgeWeight, ew);
        });
    }

    if (G.isDirected()) {
        computeReachableNodesDirected();
        G.parallelForNodes([&](node u) { margGain[u] = harmonicClosenessUBDirected(u); });
    } else {
        pred.resize(n);
        computeReachableNodesUndirected();
        G.parallelForNodes([&](node u) { margGain[u] = harmonicClosenessUBUndirected(u); });
    }
}

template <class WeightType>
GroupHarmonicClosenessImpl<WeightType>::~GroupHarmonicClosenessImpl() = default;

template <class WeightType>
void GroupHarmonicClosenessImpl<WeightType>::run() {
    group.push_back(findTopHarmonicCloseness());
#ifdef NETWORKIT_SANITY_CHECKS
    checkDistFromGroup();
#endif

    while (group.size() < k) {
        group.push_back(findNodeWithHighestMargGain());
#ifdef NETWORKIT_SANITY_CHECKS
        checkDistFromGroup();
#endif
        margGain[group.back()] = 0;
    }

    hasRun = true;
}

template <>
double GroupHarmonicClosenessImpl<count>::harmonicClosenessUBUndirected(node u) const noexcept {
    const count reachableFromU = reachableNodesInComponent[component.subsetOf(u)];
    const count degU = G->degree(u);
    double result = static_cast<double>(std::min(degU, reachableFromU));
    if (reachableFromU > degU + 1)
        result += static_cast<double>(reachableFromU - degU - 1) / 2.;
    return result;
}

template <>
double
GroupHarmonicClosenessImpl<edgeweight>::harmonicClosenessUBUndirected(node u) const noexcept {
    const count reachableFromU = reachableNodesInComponent[component.subsetOf(u)];
    if (reachableFromU <= 1)
        return 0;

    edgeweight smallestWeight = infDist;
    G->forNeighborsOf(u, [&](node v, edgeweight ew) {
        if (ew < std::min(distFromGroup[v], smallestWeight))
            smallestWeight = ew;
    });
    return 1. / smallestWeight
           + static_cast<double>(reachableFromU - 2) / (smallestWeight + minEdgeWeight);
}

template <>
double GroupHarmonicClosenessImpl<count>::harmonicClosenessUBDirected(node u) const noexcept {
    const double degU = static_cast<double>(G->degree(u));
    return degU + (static_cast<double>(reachableNodesUB[u]) - degU - 1.) / 2.;
}

template <>
double GroupHarmonicClosenessImpl<edgeweight>::harmonicClosenessUBDirected(node u) const noexcept {
    const count reachableFromU = reachableNodesUB[u];
    if (reachableFromU <= 1)
        return 0;

    const edgeweight smallestWeight =
        (*std::min_element(
             G->weightNeighborRange(u).begin(), G->weightNeighborRange(u).end(),
             [](const auto &e1, const auto &e2) -> bool { return e1.second < e2.second; }))
            .second;

    return static_cast<double>(G->degree(u)) / smallestWeight
           + static_cast<double>(reachableFromU - G->degree(u) - 1)
                 / (smallestWeight + minEdgeWeight);
}

template <class WeightType>
void GroupHarmonicClosenessImpl<WeightType>::computeReachableNodesUndirected() {
    ConnectedComponents cc(*G);
    cc.run();
    component = cc.getPartition();

    reachableNodesInComponent.clear();
    reachableNodesInComponent.reserve(cc.numberOfComponents());
    for (const auto &compSize : cc.getComponentSizes())
        reachableNodesInComponent.push_back(compSize.second);
}

template <class WeightType>
void GroupHarmonicClosenessImpl<WeightType>::computeReachableNodesDirected() {
    reachableNodesUB.resize(G->upperNodeIdBound());
    ReachableNodes rn(*G, false);
    rn.run();

    G->parallelForNodes([&reachableNodesUB = reachableNodesUB, rn](node u) {
        reachableNodesUB[u] = rn.numberOfReachableNodesUB(u);
    });
}

template <>
SSSPResult GroupHarmonicClosenessImpl<count>::prunedSSSPEmptyGroup(node source, double bestScore) {
    auto &visited = visitedGlobal[omp_get_thread_num()];
    std::fill(visited.begin(), visited.end(), false);
    visited[source] = true;

    auto &curDist = curDistGlobal[omp_get_thread_num()];
    std::fill(curDist.begin(), curDist.end(), infDist);
    curDist[source] = 0;

    const count reachableFromSource = G->isDirected()
                                          ? reachableNodesUB[source]
                                          : reachableNodesInComponent[component.subsetOf(source)];

    double curScore = 0, scoreUB = margGain[source];
    count level = 1, visitedNodes = 1;
    const count undirected = !G->isDirected();

    std::queue<node> q1, q2;
    q1.push(source);

    do {
        count nodesAtNextLevelUB = 0, prevAtNextLevelUB = 0;
        do {
            const node u = q1.front();
            q1.pop();

            for (const node w : G->neighborRange(u)) {
                if (!visited[w]) {
                    visited[w] = true;
                    curDist[w] = curDist[u] + 1;
                    q2.push(w);
                    ++visitedNodes;
                    curScore += 1. / static_cast<double>(level);
                    nodesAtNextLevelUB += G->degree(w) - undirected;
                    if (undirected)
                        pred[w] = u;
                } else if (!G->isDirected() && prevAtNextLevelUB > 0 && u != pred[w]) {
                    --prevAtNextLevelUB;
                    assert(u != source);
                    scoreUB -=
                        1. / static_cast<double>(level) - 1. / static_cast<double>(level + 1);
                    assert(scoreUB > 0);
                    if (scoreUB <= bestScore)
                        return {false, scoreUB};
                }
            }
        } while (!q1.empty());

        assert(visitedNodes <= reachableFromSource);
        nodesAtNextLevelUB = std::min(nodesAtNextLevelUB, reachableFromSource - visitedNodes);
        scoreUB = curScore
                  + static_cast<double>(nodesAtNextLevelUB) / static_cast<double>(level + 1)
                  + static_cast<double>(reachableFromSource - visitedNodes - nodesAtNextLevelUB)
                        / static_cast<double>(level + 2);

        if (scoreUB <= bestScore)
            return {false, scoreUB};
        ++level;
        std::swap(q1, q2);
    } while (!q1.empty());

    return {true, curScore};
}

template <>
SSSPResult GroupHarmonicClosenessImpl<edgeweight>::prunedSSSPEmptyGroup(node source,
                                                                        double bestScore) {
    auto &visited = visitedGlobal[omp_get_thread_num()];
    std::fill(visited.begin(), visited.end(), false);
    visited[source] = true;

    auto &curDist = curDistGlobal[omp_get_thread_num()];
    std::fill(curDist.begin(), curDist.end(), infDist);
    curDist[source] = 0;

    auto &prioQ = dijkstraHeaps[omp_get_thread_num()];

    const count reachableFromSource = G->isDirected()
                                          ? reachableNodesUB[source]
                                          : reachableNodesInComponent[component.subsetOf(source)];
    double curScore = 0, scoreUB = 0;
    count visitedNodes = 1;

    const auto exploreNeighbors = [&](node u) -> void {
        node w;
        edgeweight weight;
        for (const auto neighWeight : G->weightNeighborRange(u)) {
            std::tie(w, weight) = neighWeight;
            const edgeweight newDist = curDist[u] + weight;
            if (!visited[w]) {
                visited[w] = true;
                curDist[w] = newDist;
                prioQ.push(w);
            } else if (newDist < curDist[w]) {
                curDist[w] = newDist;
                prioQ.update(w);
            }
        }
    };

    // Explore source now to avoid "if (u != source) curScore += 1. / curDist[u];" in main
    // loop
    prioQ.clear();
    exploreNeighbors(source);

    do {
        const node u = prioQ.extract_top();
        assert(u != source);
        assert(curDist[u] > 0);
        curScore += 1. / curDist[u];
        ++visitedNodes;
        scoreUB = curScore
                  + static_cast<double>(reachableFromSource - visitedNodes)
                        / (curDist[u] + minEdgeWeight);
        if (scoreUB <= bestScore)
            return {false, scoreUB};
        exploreNeighbors(u);
    } while (!prioQ.empty());

    return {true, curScore};
}

template <>
SSSPResult GroupHarmonicClosenessImpl<count>::prunedSSSP(node source, double bestScore) {
    auto &visited = visitedGlobal[omp_get_thread_num()];
    std::fill(visited.begin(), visited.end(), false);
    visited[source] = true;

    auto &curDist = curDistGlobal[omp_get_thread_num()];
    std::fill(curDist.begin(), curDist.end(), infDist);
    curDist[source] = 0;

    auto &curNearerNodes = curNearerNodesGlobal[omp_get_thread_num()];
    curNearerNodes.clear();

    double curScore = 0, scoreUB = margGain[source];
    if (distFromGroup[source] != infDist) {
        curScore = -1. / static_cast<double>(distFromGroup[source]);
        scoreUB -= 1. / static_cast<double>(distFromGroup[source]);
    }

    count level = 1, visitedNodes = static_cast<count>(distFromGroup[source] != 1);
    const count undirected = !G->isDirected();
    const count reachableFromSource = G->isDirected()
                                          ? reachableNodesUB[source]
                                          : reachableNodesInComponent[component.subsetOf(source)];

    std::queue<node> q1, q2;
    q1.push(source);
    do {
        count nodesAtNextLevelUB = 0, prevAtNextLevelUB = 0;
        do {
            const node u = q1.front();
            q1.pop();
            curNearerNodes.push_back(u);

            for (const node w : G->neighborRange(u)) {
                if (!visited[w]) {
                    visited[w] = true;
                    ++visitedNodes;
                    if (undirected)
                        pred[w] = u;
                    if (distFromGroup[w] > level) {
                        curDist[w] = level;
                        q2.push(w);
                        nodesAtNextLevelUB += G->degree(w) - undirected;
                        curScore += 1. / static_cast<double>(level);
                        if (distFromGroup[w] != infDist)
                            curScore -= 1. / static_cast<double>(distFromGroup[w]);
                    }
                } else if (undirected && prevAtNextLevelUB && u != pred[w]) {
                    --prevAtNextLevelUB;
                    scoreUB -=
                        1. / static_cast<double>(level) - 1. / static_cast<double>(level + 1);
                    assert(scoreUB > 0);
                    if (scoreUB <= bestScore)
                        return {false, scoreUB};
                }
            }
        } while (!q1.empty());

        scoreUB =
            curScore + static_cast<double>(nodesAtNextLevelUB) / static_cast<double>(level + 1);
        if (reachableFromSource > visitedNodes + nodesAtNextLevelUB)
            scoreUB += static_cast<double>(reachableFromSource - visitedNodes - nodesAtNextLevelUB)
                       / static_cast<double>(level + 2);

        if (scoreUB <= bestScore)
            return {false, scoreUB};

        ++level;
        std::swap(q1, q2);
    } while (!q1.empty());

    return {true, curScore};
}

template <>
SSSPResult GroupHarmonicClosenessImpl<edgeweight>::prunedSSSP(node source, double bestScore) {
    auto &visited = visitedGlobal[omp_get_thread_num()];
    std::fill(visited.begin(), visited.end(), false);
    visited[source] = true;

    auto &curDist = curDistGlobal[omp_get_thread_num()];
    std::fill(curDist.begin(), curDist.end(), infDist);
    curDist[source] = 0.;

    auto &prioQ = dijkstraHeaps[omp_get_thread_num()];
    prioQ.clear();

    auto &curNearerNodes = curNearerNodesGlobal[omp_get_thread_num()];
    curNearerNodes.clear();
    curNearerNodes.push_back(source);

    const count reachableFromSource = G->isDirected()
                                          ? reachableNodesUB[source]
                                          : reachableNodesInComponent[component.subsetOf(source)];
    double curScore = 0;
    assert(distFromGroup[source] > 0);
    if (distFromGroup[source] != infDist)
        curScore = -1. / distFromGroup[source];

    count visitedNodes = 1;

    const auto exploreNeighbors = [&](node u) -> void {
        node w;
        edgeweight weight;
        for (const auto neighWeight : G->weightNeighborRange(u)) {
            std::tie(w, weight) = neighWeight;
            const edgeweight newDist = curDist[u] + weight;
            if (!visited[w]) {
                visited[w] = true;
                if (newDist < distFromGroup[w]) {
                    curDist[w] = newDist;
                    prioQ.push(w);
                }
            } else if (newDist < std::min(distFromGroup[w], curDist[w])) {
                curDist[w] = newDist;
                prioQ.update(w);
            }
        }
    };

    // Visit source now to avoid 'if (u != source)' in main loop
    exploreNeighbors(source);

    while (!prioQ.empty()) {
        const node u = prioQ.extract_top();
        ++visitedNodes;
        curNearerNodes.push_back(u);
        assert(curDist[u] < distFromGroup[u]);
        curScore += 1. / curDist[u];
        if (distFromGroup[u] != infDist)
            curScore -= 1. / distFromGroup[u];

        const double scoreUB =
            curScore + static_cast<double>(reachableFromSource - visitedNodes) / curDist[u];

        if (scoreUB <= bestScore)
            return {false, scoreUB};
        exploreNeighbors(u);
    }

    return {true, curScore};
}

template <class WeightType>
node GroupHarmonicClosenessImpl<WeightType>::findTopHarmonicCloseness() {
    node bestNode = *G->nodeRange().begin();
    double bestScore = 0;
    candidateNodesPQ.build_heap(G->nodeRange().begin(), G->nodeRange().end());

    std::atomic<bool> stop{false};

#pragma omp parallel
    {
        while (!stop.load(std::memory_order_relaxed)) {
            node u = none;

#pragma omp critical
            {
                if (candidateNodesPQ.empty()) {
                    stop.store(true, std::memory_order_relaxed);
                } else {
                    u = candidateNodesPQ.extract_top();
                    if (margGain[u] <= bestScore) {
                        stop.store(true, std::memory_order_relaxed);
                        u = none;
                    }
                }
            }

            if (u == none)
                break;

            const auto ssspResult = prunedSSSPEmptyGroup(u, bestScore);
            margGain[u] = ssspResult.score;
#pragma omp critical
            {
                if (ssspResult.complete && ssspResult.score > bestScore) {
                    bestNode = u;
                    bestScore = ssspResult.score;
                    std::swap(curDistGlobal[omp_get_thread_num()], distFromGroup);
#ifdef NETWORKIT_SANITY_CHECKS
                    group.push_back(bestNode);
                    checkDistFromGroup();
                    group.clear();
#endif // NETWORKIT_SANITY_CHECKS
                }
            }
        }
    }

    if (!G->isDirected() && !G->isWeighted())
        reachableNodesInComponent[component.subsetOf(bestNode)] -= G->degree(bestNode) + 1;
    return bestNode;
}

template <class WeightType>
node GroupHarmonicClosenessImpl<WeightType>::findNodeWithHighestMargGain() {
    node bestNode = none;
    double bestScore = -std::numeric_limits<double>::max();
    if (!G->isDirected())
        G->forNodes([&](node u) {
            if (distFromGroup[u] > 0)
                margGain[u] = std::min(margGain[u], harmonicClosenessUBUndirected(u));
        });

    candidateNodesPQ.build_heap(G->nodeRange().begin(), G->nodeRange().end());
    for (node u : group)
        candidateNodesPQ.remove(u);

    std::atomic<bool> stop{false};
#pragma omp parallel
    {
        while (!stop.load(std::memory_order_relaxed)) {
            node u = none;
#pragma omp critical
            {
                if (candidateNodesPQ.empty()) {
                    stop.store(true, std::memory_order_relaxed);
                } else {
                    u = candidateNodesPQ.extract_top();
                    if (margGain[u] <= bestScore) {
                        stop.store(true, std::memory_order_relaxed);
                        u = none;
                    }
                }
            }

            if (u == none)
                break;

            const auto ssspResult = prunedSSSP(u, bestScore);
            margGain[u] = ssspResult.score;

#pragma omp critical
            {
                if (ssspResult.complete && ssspResult.score > bestScore) {
                    bestNode = u;
                    bestScore = ssspResult.score;
                    std::swap(curDistGlobal[omp_get_thread_num()], curBestDist);
                    std::swap(curNearerNodesGlobal[omp_get_thread_num()], nearerNodes);
                }
            }
        }
    }

    assert(bestNode != none);
    for (const node u : nearerNodes) {
        assert(distFromGroup[u] > curBestDist[u]);
        distFromGroup[u] = curBestDist[u];
        if (!G->isWeighted() && !G->isDirected() && distFromGroup[u] > 1 && curBestDist[u] <= 1) {
            assert(reachableNodesInComponent[component.subsetOf(u)]);
            --reachableNodesInComponent[component.subsetOf(u)];
        }
    }

    return bestNode;
}

template <>
double GroupHarmonicClosenessImpl<count>::scoreOfGroup(const Graph &graph,
                                                       const std::vector<node> &inputGroup) {
    double score = 0.;
    Traversal::BFSfrom(graph, inputGroup.begin(), inputGroup.end(),
                       [&score](node, edgeweight dist) {
                           if (dist > 0.)
                               score += 1. / dist;
                       });
    return score;
}

template <>
double GroupHarmonicClosenessImpl<edgeweight>::scoreOfGroup(const Graph &graph,
                                                            const std::vector<node> &inputGroup) {
    double score = 0.;
    Traversal::DijkstraFrom(graph, inputGroup.begin(), inputGroup.end(),
                            [&score](node, edgeweight dist) {
                                if (dist > 0.)
                                    score += 1. / dist;
                            });
    return score;
}

#ifdef NETWORKIT_SANITY_CHECKS
template <>
void GroupHarmonicClosenessImpl<count>::checkDistFromGroup() const {
    Traversal::BFSfrom(
        *G, group.begin(), group.end(),
        [&distFromGroup = distFromGroup](node x, count d) { assert(distFromGroup[x] == d); });
}

template <>
void GroupHarmonicClosenessImpl<edgeweight>::checkDistFromGroup() const {
    Traversal::DijkstraFrom(
        *G, group.begin(), group.end(),
        [&distFromGroup = distFromGroup](node x, edgeweight d) { assert(distFromGroup[x] == d); });
}
#endif // NETWORKIT_SANITY_CHECKS

} // namespace

using GroupHarmonicClosenessUnweighted = GroupHarmonicClosenessImpl<count>;
using GroupHarmonicClosenessWeighted = GroupHarmonicClosenessImpl<edgeweight>;

GroupHarmonicCloseness::GroupHarmonicCloseness(const Graph &G, count k) : weighted(G.isWeighted()) {
    if (weighted)
        impl = std::make_unique<GroupHarmonicClosenessWeighted>(G, k);
    else
        impl = std::make_unique<GroupHarmonicClosenessUnweighted>(G, k);
}

void GroupHarmonicCloseness::run() {
    impl->run();
    hasRun = true;
}

const std::vector<node> &GroupHarmonicCloseness::groupMaxHarmonicCloseness() const {
    assureFinished();
    return impl->group;
}

double GroupHarmonicCloseness::scoreOfGroup(const Graph &graph,
                                            const std::vector<node> &inputGroup) {
    if (graph.isWeighted())
        return GroupHarmonicClosenessWeighted::scoreOfGroup(graph, inputGroup);
    else
        return GroupHarmonicClosenessUnweighted::scoreOfGroup(graph, inputGroup);
}

} // namespace NetworKit
