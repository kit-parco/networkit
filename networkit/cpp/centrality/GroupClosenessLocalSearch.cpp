
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <queue>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <utility>

#ifdef __AVX2__
#include <immintrin.h>
#include <networkit/auxiliary/AlignedAllocator.hpp>
#endif // __AVX2__

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/VectorComparator.hpp>
#include <networkit/centrality/GroupClosenessGrowShrink.hpp>
#include <networkit/centrality/GroupClosenessLocalSearch.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/Dijkstra.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>
#include <tlx/container/d_ary_heap.hpp>
#include <tlx/unused.hpp>

namespace NetworKit {

namespace {

template <class WeightType>
class GroupClosenessLocalSearchImpl final
    : public GroupClosenessLocalSearch::GroupClosenessLocalSearchInterface {

    static constexpr count K = 16; // Number of packed repetitions using AVX
    static constexpr float maxInt16 = static_cast<float>(std::numeric_limits<uint16_t>::max());
    static constexpr WeightType invalidDecr = static_cast<WeightType>(-1);
    static constexpr WeightType infDistance = std::numeric_limits<WeightType>::max();

public:
    template <class InputIt>
    GroupClosenessLocalSearchImpl(const Graph &G, InputIt initGroupFirst, InputIt initGroupLast,
                                  bool runGrowShrink, count maxIterations);

    void run() override;

private:
    const Graph *G;
    const bool runGrowShrink;
    const count maxIterations;

    // Shortest and 2nd shortest distance between the group to all the nodes
    std::vector<WeightType> distance, distance2;
    // Shortest and 2nd shortest distances before a swap
    std::vector<WeightType> oldDistance, oldDistance2;
    // Distance vectors used by SSSPs to avoid dynamic memory allocations
    std::vector<std::vector<WeightType>> distanceGlobal, distance2Global;
    std::vector<WeightType> farnessIncrease;

    // First and second nearest node in the group to each other node
    std::vector<node> nearest, nearest2;
    // Same as above, but use thread-local vector when running in parallel
    std::vector<std::vector<node>> nearestGlobal, nearest2Global;
    // First and second nearest node in the group to each other node before a swap
    std::vector<node> oldNearest, oldNearest2;

    // Used for graph exploration to avoid dynamic memory allocations
    std::vector<std::vector<bool>> visitedGlobal;

    // Stores at index 'i' the number of nodes at distance 'i' from the group
    std::vector<node> nodesAtDistance;
    // All nodes sorted by their distance from the group
    std::vector<node> nodesSortedByDistance;
    // At each index 'i', farness due to the vertices farther than the vertex at distance
    // nodesSortedByDistance[i]
    std::vector<WeightType> sumOfDistancesGreaterThan;
    // When doing a pruned SSSP to evaluate a swap, stores at index 'i' how many unexplored vertices
    // are at distance 'i' from the source (one vector per thread).
    std::vector<std::vector<node>> nodesLeftToExploreAtDistanceGlobal;
    // Each thread independently evaluates the decrease in farness after the insertion of a vertex
    // to the group. For each thread we store a pair <u, decrement>, i.e., the decrement in farness
    // due to the insertion of vertex 'u' to the group.
    std::vector<std::pair<node, WeightType>> decreasePerThread;

    // Stack used to compute the DAG with a SSSP from the group
    std::vector<node> dagStack;
    // Nodes to be added to the group sorted by increasing estimated decrease in farness
    std::vector<node> candidatesToAddVec;
    // Estimate decrease in farness for each candidate node
    std::vector<float> approxFarnessDecrease;
    // Used by the randomized algorithm that estimates the reachability set size of a vertex in a
    // DAG
    std::vector<count> sumOfMins;

    count diam, stackSize;

    WeightType groupFarness;

    tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<WeightType>> heap, heap2;
    std::vector<tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<WeightType>>> heaps;
    tlx::d_ary_heap<node, 2, Aux::GreaterInVector<float>> candidatesToAdd;
    tlx::d_ary_heap<node, 2, Aux::LessInVector<WeightType>> candidatesToRemove;

    std::vector<std::reference_wrapper<std::mt19937_64>> urngs;
#ifdef __AVX2__
    union RandItem {
        uint16_t items[K];
        __m256i vec;
    };
    std::vector<RandItem, AlignedAllocator<RandItem, sizeof(RandItem)>> randVec;
#else
    std::vector<uint16_t> randVec;
#endif // __AVX2__

    std::vector<std::uniform_int_distribution<uint32_t>> intDistributions;

    void init();
    count estimateDiamDirected();
    void computeDistances();
    void traversalFromGroup();
    tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<WeightType>> &getPQ();
    void dijkstra();
    bool findAndSwap();
    void computeFarnessIncrease();
    void updateAfterRemoval(node u, bool updateNodesAtDistance = false);
    node estimateFarnessDecrease(bool getTopOnly = false);
    WeightType computeFarnessDecrease(node v, WeightType increaseU = 0, bool pruned = false,
                                      bool parallel = false);
    void computeDAG();
    void initRandomVec();

    // Stopping condition of the algorithm: determines whether the decrease in farness achieved
    // is enough to proceed. If not, the algorithm stops.
    bool insuffcientDecrease(WeightType decrease, WeightType increase) const noexcept {
        const count k = group.size(), n = G->numberOfNodes();
        return decrease <= increase
               || static_cast<double>(decrease - increase) / static_cast<double>(groupFarness)
                      < 1. / static_cast<double>(k * (n - k));
    }

#ifdef NETWORKIT_SANITY_CHECKS
    void runSanityCheck(node v = none, bool parallel = false) const;
#endif // NETWORKIT_SANITY_CHECKS
};

template <class WeightType>
template <class InputIt>
GroupClosenessLocalSearchImpl<WeightType>::GroupClosenessLocalSearchImpl(const Graph &G,
                                                                         InputIt initGroupFirst,
                                                                         InputIt initGroupLast,
                                                                         bool runGrowShrink,
                                                                         count maxIterations)
    : GroupClosenessLocalSearchInterface(initGroupFirst, initGroupLast), G(&G),
      runGrowShrink(runGrowShrink), maxIterations(maxIterations),
      heap(Aux::LessInVector<WeightType>(distance)),
      heap2(Aux::LessInVector<WeightType>(distance2)),
      candidatesToAdd(Aux::GreaterInVector<float>(approxFarnessDecrease)),
      candidatesToRemove(Aux::LessInVector<WeightType>(farnessIncrease)) {

    if (group.empty())
        throw std::runtime_error("Error, empty group.");
}

template <class WeightType>
count GroupClosenessLocalSearchImpl<WeightType>::estimateDiamDirected() {
    auto &visited = visitedGlobal[omp_get_thread_num()];

    count estDiam = 0;

    // Farthest node from u in #of hops (weights are ignored)
    const auto farthestFrom = [&](node u) -> node {
        std::fill(visited.begin(), visited.end(), false);
        visited[u] = true;

        std::queue<node> q1, q2;
        q1.push(u);
        node farthest = u;

        node dist = 1;
        do {
            do {
                const node v = q1.front();
                q1.pop();
                G->forNeighborsOf(v, [&](const node w) {
                    if (!visited[w]) {
                        visited[w] = true;
                        farthest = u;
                        q2.push(w);
                    }
                });

            } while (!q1.empty());

            estDiam = std::max(estDiam, dist);
            ++dist;
            std::swap(q1, q2);
        } while (!q1.empty());

        return farthest;
    };

    node source =
        *std::max_element(G->nodeRange().begin(), G->nodeRange().end(),
                          [G = G](node u, node v) { return G->degree(u) < G->degree(v); });
    // Preliminary experiments show that 10 BFSs are enough to get a good estimate of the diameter
    for (int i = 0; i < 10; ++i)
        source = farthestFrom(source);

    return estDiam;
}

template <class WeightType>
void GroupClosenessLocalSearchImpl<WeightType>::init() {
    const auto n = G->upperNodeIdBound();

    distance.assign(n, 0);
    distance2.assign(n, 0);
    distanceGlobal.assign(omp_get_max_threads(), std::vector<WeightType>(n));
    distance2Global.assign(omp_get_max_threads(), std::vector<WeightType>(n));

    visitedGlobal.assign(omp_get_max_threads(), std::vector<bool>(n));

    nearest.assign(n, none);
    nearest2.assign(n, none);
    oldNearest.assign(n, none);
    oldNearest2.assign(n, none);
    nearestGlobal.assign(omp_get_max_threads(), std::vector<node>(n, none));
    nearest2Global.assign(omp_get_max_threads(), std::vector<node>(n, none));

    decreasePerThread.resize(omp_get_max_threads(), std::make_pair(none, 0));

    approxFarnessDecrease.assign(n, 0);
    farnessIncrease.assign(n, 0);

    heap.reserve(n);
    candidatesToAdd.reserve(n);
    candidatesToRemove.reserve(n);

    candidatesToAddVec.reserve(n - group.size());
    dagStack.assign(n, 0);
    sumOfMins.assign(n, 0);

    if (G->isDirected()) {
        diam = estimateDiamDirected();
    } else {
        Diameter diamAlgo(*G, DiameterAlgo::estimatedRange, 0.1);
        diamAlgo.run();
        diam = diamAlgo.getDiameter().second;
    }

    if (G->isWeighted()) {
        nodesSortedByDistance.assign(G->nodeRange().begin(), G->nodeRange().end());
        sumOfDistancesGreaterThan.resize(n);
        heap2.reserve(n);
        heaps.reserve(omp_get_max_threads());
        for (omp_index i = 0; i < static_cast<omp_index>(omp_get_max_threads()); ++i) {
            heaps.emplace_back(distanceGlobal[i]);
            heaps.back().reserve(n);
        }
    } else {
        nodesAtDistance.resize(diam + 1);
        nodesLeftToExploreAtDistanceGlobal.resize(omp_get_max_threads(),
                                                  std::vector<node>(diam + 1));
    }

    urngs.reserve(omp_get_max_threads());
    for (omp_index i = 0; i < static_cast<omp_index>(omp_get_max_threads()); ++i)
        urngs.emplace_back(Aux::Random::getURNG());

#ifdef __AVX2__
    randVec.resize(n);
#else
    randVec.resize(K * n);
#endif // __AVX2__

    intDistributions.resize(omp_get_max_threads());
}

template <>
tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<count>> &
GroupClosenessLocalSearchImpl<count>::getPQ() {
    return heap;
}

template <>
tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<edgeweight>> &
GroupClosenessLocalSearchImpl<edgeweight>::getPQ() {
    return heap2;
}

template <class WeightType>
void GroupClosenessLocalSearchImpl<WeightType>::initRandomVec() {
    G->parallelForNodes([&](node u) {
        // Avoid to generate numbers for nodes in the group
        if (distance[u] > 0) {
            auto tid = omp_get_thread_num();
            auto &curUrng = urngs[tid].get();
            auto &distr = intDistributions[tid];
#ifdef __AVX2__
            // Generating two 16-bit random integers per time
            for (index j = 0; j < K; j += 2) {
                const auto x = distr(curUrng);
                randVec[u].items[j] = static_cast<uint16_t>(x);
                randVec[u].items[j + 1] = static_cast<uint16_t>(x >> K);
            }
            randVec[u].vec = *(__m256i *)(&randVec[u].items[0]);
#else
            // Generating two 16-bit random integers per time
            for (index j = 0; j < K; j += 2) {
                const auto x = distr(curUrng);
                randVec[K * u + j] = static_cast<uint16_t>(x);
                randVec[K * u + j + 1] = static_cast<uint16_t>(x >> K);
            }
#endif // __AVX2__
        }
    });
}

// Unweighted case: perform a BFS
template <>
void GroupClosenessLocalSearchImpl<count>::traversalFromGroup() {
    groupFarness = 0;
    auto &visited = visitedGlobal[omp_get_thread_num()];
    std::fill(visited.begin(), visited.end(), false);

    std::queue<node> q;

    for (node u : group) {
        q.push(u);
        distance[u] = 0;
        visited[u] = true;
        nearest[u] = u;
    }

    do {
        const node x = q.front();
        q.pop();
        groupFarness += distance[x];

        G->forNeighborsOf(x, [&](node y) {
            if (!visited[y]) {
                distance[y] = distance[x] + 1;
                nearest[y] = nearest[x];
                visited[y] = true;
                q.push(y);
            }
        });

    } while (!q.empty());
}

// Weighted case: Dijkstra
template <>
void GroupClosenessLocalSearchImpl<edgeweight>::traversalFromGroup() {
    groupFarness = 0;
    assert(heap.empty());
    auto &visited = visitedGlobal[omp_get_thread_num()];
    std::fill(visited.begin(), visited.end(), false);

    for (node u : group) {
        heap.push(u);
        distance[u] = 0;
        visited[u] = true;
        nearest[u] = u;
    }

    do {
        const node x = heap.extract_top();
        groupFarness += distance[x];

        G->forNeighborsOf(x, [&](node y, edgeweight weight) {
            if (!visited[y] || distance[y] > distance[x] + weight) {
                distance[y] = distance[x] + weight;
                nearest[y] = nearest[x];
                visited[y] = true;
                heap.update(y);
            }
        });

    } while (!heap.empty());
}

template <class WeightType>
void GroupClosenessLocalSearchImpl<WeightType>::dijkstra() {
    auto &pq = getPQ();
    assert(!pq.empty());

    auto &visited = visitedGlobal[omp_get_thread_num()];

    do {
        const node x = pq.extract_top();

        G->forNeighborsOf(x, [&](node y, edgeweight weight) {
            if (!visited[y] || distance2[y] > distance2[x] + weight) {
                distance2[y] = distance2[x] + weight;
                nearest2[y] = nearest2[x];
                pq.update(y);
                visited[y] = true;
            }
        });

    } while (!pq.empty());
}

template <class WeightType>
void GroupClosenessLocalSearchImpl<WeightType>::computeDistances() {
    // Compute 1st shortest distances
    traversalFromGroup();

    // Compute 2nd shortest distances
    // Step 1: find 'boundary' nodes (sources of Dijkstra)
    auto &visited = visitedGlobal[omp_get_thread_num()];
    std::fill(visited.begin(), visited.end(), false);
    auto &pq = getPQ();

    const auto processEdge = [&](node x, node y, WeightType weight) -> void {
        if (!visited[y] || distance2[y] > distance[x] + weight) {
            distance2[y] = distance[x] + weight;
            visited[y] = true;
            nearest2[y] = nearest[x];
            pq.update(y);
        }
    };

    G->forEdges([&](node x, node y, edgeweight weight) {
        if (nearest[x] != nearest[y]) {
            processEdge(x, y, static_cast<WeightType>(weight));
            if (!G->isDirected())
                processEdge(y, x, static_cast<WeightType>(weight));
        }
    });

    // Step 2: Explore rest of the graph with Dijkstra
    dijkstra();
}

template <class WeightType>
void GroupClosenessLocalSearchImpl<WeightType>::computeFarnessIncrease() {
    for (node u : group)
        farnessIncrease[u] = 0;

    G->forNodes([&](node u) {
        assert(group.find(nearest[u]) != group.end());
        farnessIncrease[nearest[u]] += distance2[u] - distance[u];
    });
}

// Unweighted graphs: use a BFS
template <>
void GroupClosenessLocalSearchImpl<count>::computeDAG() {
    auto &visited = visitedGlobal[omp_get_thread_num()];
    std::fill(visited.begin(), visited.end(), false);

    std::queue<node> q;
    for (node u : group) {
        q.push(u);
        visited[u] = true;
    }

    // Recompute the DAG
    do {
        const node x = q.front();
        q.pop();

        bool isLeaf = true;
        G->forNeighborsOf(x, [&](node y) {
            if (!visited[y]) {
                visited[y] = true;
                isLeaf = false;
                q.push(y);
            }
        });

        if (distance[x] > 0 && (!isLeaf || distance[x] == 1))
            dagStack[stackSize++] = x;

    } while (!q.empty());
}

// Weighted graphs: use Dijkstra
template <>
void GroupClosenessLocalSearchImpl<edgeweight>::computeDAG() {
    assert(heap.empty());
    auto &visited = visitedGlobal[omp_get_thread_num()];
    std::fill(visited.begin(), visited.end(), false);

    for (node u : group) {
        heap.push(u);
        visited[u] = true;
    }

    // Recompute the DAG
    do {
        const node x = heap.extract_top();

        G->forNeighborsOf(x, [&](node y, edgeweight weight) {
            if (!visited[y] || (distance[y] > distance[x] + weight)) {
                visited[y] = true;
                heap.update(y);
            }
        });

        if (distance[x] > 0)
            dagStack[stackSize++] = x;

    } while (!heap.empty());
}

template <class WeightType>
node GroupClosenessLocalSearchImpl<WeightType>::estimateFarnessDecrease(bool getTopOnly) {
    stackSize = 0;
    computeDAG();

    // Generating all the necessary random numbers
    initRandomVec();

    // Do 16 packed repetitions at once
    for (count i = 0; i < stackSize; ++i) {
        const node x = dagStack[stackSize - 1 - i];

#ifdef __AVX2__
        // 16 randomly generated integers;
        __m256i &x1 = randVec[x].vec;
        // Pulling leaves
        G->forNeighborsOf(x, [&](node y, edgeweight weight) {
            if (distance[y] == distance[x] + weight) {
                const __m256i &y1 = randVec[y].vec;
                x1 = _mm256_min_epu16(x1, y1);
            }
        });
        *(__m256i *)(&randVec[x].items) = x1;
#else
        G->forNeighborsOf(x, [&](node y, edgeweight weight) {
            if (distance[y] == distance[x] + weight)
                for (index i = 0; i < K; ++i)
                    randVec[K * x + i] = std::min(randVec[K * x + i], randVec[K * y + i]);
        });

#endif // __AVX2__

        if (distance[x] > 0) {
            sumOfMins[x] = 0;
            for (index j = 0; j < K; ++j)
#ifdef __AVX2__
                sumOfMins[x] += randVec[x].items[j];
#else
                sumOfMins[x] = randVec[K * x + j];
#endif // __AVX2__

            if (sumOfMins[x] == 0)
                sumOfMins[x] = 1;
        }
    }

    if (!getTopOnly)
        candidatesToAddVec.clear();

    node v = none;
    float bestEstimate = -1.f;
    G->forNodes([&](const node x) {
        if (distance[x] > 0 && sumOfMins[x] > 0) {
            const float estimate = static_cast<float>(distance[x])
                                   * (16.f / (static_cast<float>(sumOfMins[x]) / maxInt16) - 1.f);
            if (getTopOnly) {
                if (estimate > bestEstimate) {
                    bestEstimate = estimate;
                    v = x;
                }
            } else if (G->degree(x) > 1) {
                approxFarnessDecrease[x] = estimate;
                candidatesToAddVec.push_back(x);
            }
        } else if (!getTopOnly)
            approxFarnessDecrease[x] = -1.f;
    });

    return v;
}

// Unweighted graphs
template <>
count GroupClosenessLocalSearchImpl<count>::computeFarnessDecrease(node v, count increaseU,
                                                                   bool pruned, bool parallel) {
    auto &curDistance = parallel ? distanceGlobal[omp_get_thread_num()] : distance;
    auto &curDistance2 = parallel ? distance2Global[omp_get_thread_num()] : distance2;
    auto &curNearest = parallel ? nearestGlobal[omp_get_thread_num()] : nearest;
    auto &curNearest2 = parallel ? nearest2Global[omp_get_thread_num()] : nearest2;

    count decr = 0;
    curDistance2[v] = curDistance[v];
    curNearest2[v] = curNearest[v];
    curDistance[v] = 0;
    curNearest[v] = v;

    std::queue<node> q1, q2;
    q1.push(v);

    auto &visited = visitedGlobal[omp_get_thread_num()];
    std::fill(visited.begin(), visited.end(), false);
    visited[v] = true;

    auto &nodesLeftToExploreAtDistance = nodesLeftToExploreAtDistanceGlobal[omp_get_thread_num()];
    if (pruned) {
        std::copy(nodesAtDistance.begin(), nodesAtDistance.end(),
                  nodesLeftToExploreAtDistance.begin());
        --nodesLeftToExploreAtDistance[curDistance[v]];
    }

    count level = 0;
    do {
        count nextLevelDecr = 0;

        do {
            const node x = q1.front();
            q1.pop();

            decr += static_cast<count>(curNearest[x] == v) * (curDistance2[x] - curDistance[x]);

            G->forNeighborsOf(x, [&](const node y) {
                if (visited[y])
                    return;

                const count distY = curDistance[y];

                if (pruned)
                    --nodesLeftToExploreAtDistance[distY];

                if (distY > curDistance[x] + 1) {

                    curDistance2[y] = distY;
                    curNearest2[y] = curNearest[y];
                    curDistance[y] = curDistance[x] + 1;
                    curNearest[y] = v;
                    q2.push(y);
                    if (pruned)
                        nextLevelDecr += static_cast<count>(distY - curDistance[y]);
                } else if (curNearest[x] == v && curDistance2[y] > curDistance[x] + 1) {
                    curDistance2[y] = curDistance[x] + 1;
                    curNearest2[y] = v;
                    q2.push(y);
                } else if (curNearest2[x] == v && curDistance2[y] > curDistance2[x] + 1) {
                    curDistance2[y] = curDistance2[x] + 1;
                    curNearest2[y] = v;
                    q2.push(y);
                }

                visited[y] = true;
            });

        } while (!q1.empty());

        if (pruned) {
            // Upper bound to decrease:
            // decrease due to visited vertices +
            // decrease due to enqueued vertices (in the next level of the BFS) +
            // (vertices at distance level + 2 do not carry any decrease because they stay at the
            // same distance) +
            // vertices farther than current level + 2 are assumed to be reachable from v in
            // level + 2 hops, so they could bring a decrease of (dist(x, S) - level - 2).
            // e.g., if nodes at distance (level + 3) can be reached in (level + 2) hops,
            // will bring a decrease of 1
            count decrUB = decr + nextLevelDecr;
            for (count dist = level + 3; dist < diam; ++dist)
                decrUB += nodesLeftToExploreAtDistance[dist] * (dist - level - 2);

            // Prune
            if (insuffcientDecrease(decrUB, increaseU))
                return invalidDecr;
        }

        ++level;
        std::swap(q1, q2);
    } while (!q1.empty());

    return decr;
}

// Weighted graphs
template <>
edgeweight
GroupClosenessLocalSearchImpl<edgeweight>::computeFarnessDecrease(node v, edgeweight increaseU,
                                                                  bool pruned, bool parallel) {
    auto &curDistance = parallel ? distanceGlobal[omp_get_thread_num()] : distance;
    auto &curDistance2 = parallel ? distance2Global[omp_get_thread_num()] : distance2;
    auto &curNearest = parallel ? nearestGlobal[omp_get_thread_num()] : nearest;
    auto &curNearest2 = parallel ? nearest2Global[omp_get_thread_num()] : nearest2;
    auto &curHeap = parallel ? heaps[omp_get_thread_num()] : heap;

    curDistance2[v] = curDistance[v];
    curNearest2[v] = curNearest[v];
    curDistance[v] = 0;
    curNearest[v] = v;

    curHeap.push(v);
    edgeweight decr = 0;

    do {
        const node x = curHeap.extract_top();
        decr += static_cast<edgeweight>(curNearest[x] == v) * (curDistance2[x] - curDistance[x]);

        G->forNeighborsOf(x, [&](const node y, const edgeweight weight) {
            if (curDistance[y] > curDistance[x] + weight) {
                // Nearest nodes could have been already updated
                if (curNearest[y] != v) {
                    curNearest2[y] = curNearest[y];
                    curNearest[y] = v;
                    curDistance2[y] = curDistance[y];
                }
                curDistance[y] = curDistance[x] + weight;
                curHeap.update(y);
            } else if (curNearest[x] == v && curNearest[y] != v
                       && curDistance2[y] > curDistance[x] + weight) {
                curDistance2[y] = curDistance[x] + weight;
                curNearest2[y] = v;
                curHeap.update(y);
            } else if (curNearest2[x] == v && curNearest[y] != v
                       && curDistance2[y] > curDistance2[x] + weight) {
                curDistance2[y] = curDistance2[x] + weight;
                curNearest2[y] = v;
                curHeap.update(y);
            }
        });

        if (curHeap.empty())
            break;

        if (pruned) {
            const count n = G->numberOfNodes();
            const node top = curHeap.top();
            const edgeweight topDist = distance[top];
            assert(curDistance2[top] >= topDist);
            edgeweight decrUB =
                decr
                + static_cast<edgeweight>(curNearest[top] == v) * (curDistance2[top] - topDist);
            const count nextIdx =
                (std::lower_bound(
                     nodesSortedByDistance.begin(), nodesSortedByDistance.end(), top,
                     [&](const node a, const node b) { return distance[a] < distance[b]; })
                 - nodesSortedByDistance.begin())
                + 1;

            if (nextIdx < n) {
                const auto nodesLeft = static_cast<edgeweight>(n - nextIdx);
                decrUB += sumOfDistancesGreaterThan[nextIdx] - nodesLeft * topDist;
                assert(decrUB >= decr);
            }

            if (insuffcientDecrease(decrUB, increaseU)) {
                curHeap.clear();
                return invalidDecr;
            }
        }
    } while (true);

    return decr;
}

// Unweighted graphs
template <>
void GroupClosenessLocalSearchImpl<count>::updateAfterRemoval(node u, bool updateNodesAtDistance) {
    G->forNodes([&](node x) {
        tlx::unused(x);
        assert(nearest[x] != nearest2[x]);
    });

    if (updateNodesAtDistance)
        std::fill(nodesAtDistance.begin(), nodesAtDistance.end(), 0);

    // Updating (1st) shortest distances and setting 2nd shortest distances to infinity
    G->forNodes([&](node x) {
        if (nearest[x] == u) {
            nearest[x] = nearest2[x];
            distance[x] = distance2[x];
            nearest2[x] = none;
            distance2[x] = infDistance;
            assert(nearest[x] != u);
        } else if (nearest2[x] == u) {
            nearest2[x] = none;
            distance2[x] = infDistance;
        }
        if (updateNodesAtDistance)
            ++nodesAtDistance[distance[x]];
    });

    auto &pq = getPQ();

    // Update the distance of y
    const auto processEdge = [&](node x, node y, count weight) -> void {
        if (nearest[x] != nearest[y]) {
            if (distance2[y] > distance[x] + weight) {
                distance2[y] = distance[x] + weight;
                nearest2[y] = nearest[x];
                pq.update(y);
            }
            // Checking 2nd distances
        } else if (distance2[x] < infDistance && distance2[y] > distance2[x] + weight) {
            distance2[y] = distance2[x] + weight;
            nearest2[y] = nearest2[x];
            pq.update(y);
            assert(nearest2[x] != nearest[y]);
        }
    };

    G->forEdges([&](node x, node y, edgeweight weight) {
        processEdge(x, y, static_cast<count>(weight));
        if (!G->isDirected())
            processEdge(y, x, static_cast<count>(weight));
    });

    while (!pq.empty()) {
        const node x = pq.extract_top();
        G->forNeighborsOf(x, [&](node y, edgeweight weight) {
            if (nearest[y] != nearest2[x]
                && distance2[y] > distance2[x] + static_cast<count>(weight)) {
                distance2[y] = distance2[x] + static_cast<count>(weight);
                nearest2[y] = nearest2[x];
                pq.update(y);
            }
        });
    }

    G->forNodes([&](node x) {
        tlx::unused(x);
        assert(nearest[x] != u);
        assert(nearest2[x] != u);
    });
}

// Weighted graphs
template <>
void GroupClosenessLocalSearchImpl<edgeweight>::updateAfterRemoval(node u,
                                                                   bool updateNodesAtDistance) {
    G->forNodes([&](node x) {
        tlx::unused(x);
        assert(nearest[x] != nearest2[x]);
    });

    // Updating (1st) shortest distances and setting 2nd shortest distances to infinity
    G->forNodes([&](const node x) {
        if (nearest[x] == u) {
            nearest[x] = nearest2[x];
            distance[x] = distance2[x];
            nearest2[x] = none;
            distance2[x] = infDistance;
            assert(nearest[x] != u);
        } else if (nearest2[x] == u) {
            nearest2[x] = none;
            distance2[x] = infDistance;
        }
    });

    auto &pq = getPQ();

    // Update the distance of y
    const auto process = [&](node x, node y, edgeweight weight) -> void {
        if (nearest[x] != nearest[y]) {
            if (distance2[y] > distance[x] + weight) {
                distance2[y] = distance[x] + weight;
                nearest2[y] = nearest[x];
                pq.update(y);
            }
            // Checking 2nd distances
        } else if (distance2[x] < infDistance && distance2[y] > distance2[x] + weight) {
            distance2[y] = distance2[x] + weight;
            nearest2[y] = nearest2[x];
            pq.update(y);
            assert(nearest2[x] != nearest[y]);
        }
    };

    G->forEdges([&](node x, node y, edgeweight weight) {
        process(x, y, weight);
        if (!G->isDirected())
            process(y, x, weight);
    });

    while (!pq.empty()) {
        const node x = pq.extract_top();
        G->forNeighborsOf(x, [&](node y, edgeweight weight) {
            if (nearest[y] != nearest2[x] && distance2[y] > distance2[x] + weight) {
                distance2[y] = distance2[x] + weight;
                nearest2[y] = nearest2[x];
                pq.update(y);
            }
        });
    }

    if (updateNodesAtDistance) {
        Aux::Parallel::sort(
            nodesSortedByDistance.begin(), nodesSortedByDistance.end(),
            [&distance = distance](node x, node y) { return distance[x] < distance[y]; });

        sumOfDistancesGreaterThan[G->numberOfNodes() - 1] = distance[nodesSortedByDistance.back()];
        for (node i = G->numberOfNodes() - 1; i; --i)
            sumOfDistancesGreaterThan[i - 1] =
                sumOfDistancesGreaterThan[i] + distance[nodesSortedByDistance[i - 1]];
    }
}

template <class WeightType>
bool GroupClosenessLocalSearchImpl<WeightType>::findAndSwap() {

    // Copy distances/nearest vertices from shared vectors to thread-local storage
    const auto copyStatusToThread = [&]() -> void {
        const auto thread = omp_get_thread_num();
        std::copy(distance.begin(), distance.end(), distanceGlobal[thread].begin());
        std::copy(distance2.begin(), distance2.end(), distance2Global[thread].begin());
        std::copy(nearest.begin(), nearest.end(), nearestGlobal[thread].begin());
        std::copy(nearest2.begin(), nearest2.end(), nearest2Global[thread].begin());
    };

    // Store distances/nearest vertices from thread-local storage to shared vectors
    const auto storeStatusFromThread = [&](const auto thread) -> void {
        std::copy(distanceGlobal[thread].begin(), distanceGlobal[thread].end(), distance.begin());
        std::copy(distance2Global[thread].begin(), distance2Global[thread].end(),
                  distance2.begin());
        std::copy(nearestGlobal[thread].begin(), nearestGlobal[thread].end(), nearest.begin());
        std::copy(nearest2Global[thread].begin(), nearest2Global[thread].end(), nearest2.begin());
    };

    // Update the farness increase due to the removal of every node in the group
    computeFarnessIncrease();
    candidatesToRemove.build_heap(group.begin(), group.end());

    // Iterate over all vertices in the group
    do {
        const node u = candidatesToRemove.extract_top();
        const WeightType increaseU = farnessIncrease[u];

        group.erase(u);
        updateAfterRemoval(u, true);
#ifdef NETWORKIT_SANITY_CHECKS
        runSanityCheck();
#endif // NETWORKIT_SANITY_CHECKS

        // Approximate the farness decrease for each vertex by estimating their BFS DAG
        // size
        estimateFarnessDecrease();
        candidatesToAdd.build_heap(candidatesToAddVec.begin(), candidatesToAddVec.end());

        std::atomic<bool> stop{false};
        std::atomic<int> threadToSelect{-1};

#pragma omp parallel
        {
            while (!stop.load(std::memory_order_relaxed)) {

                node v = none;
#pragma omp critical
                {
                    if (candidatesToAdd.empty()) {
                        stop.store(true, std::memory_order_relaxed);
                    } else {
                        v = candidatesToAdd.extract_top();
                        if (approxFarnessDecrease[v] <= 0.f) {
                            stop.store(true, std::memory_order_relaxed);
                            v = none;
                        }
                    }
                }

                if (v == none)
                    break;

                copyStatusToThread();
#ifdef NETWORKIT_SANITY_CHECKS
                runSanityCheck(none, true);
#endif // NETWORKIT_SANITY_CHECKS

                const WeightType farnessDecrease = computeFarnessDecrease(v, increaseU, true, true);

#ifdef NETWORKIT_SANITY_CHECKS
                if (farnessDecrease != invalidDecr)
                    runSanityCheck(v, true);
#endif // NETWORKIT_SANITY_CHECKS

                if (farnessDecrease == invalidDecr
                    || insuffcientDecrease(farnessDecrease, increaseU))
                    continue;

                stop.store(true, std::memory_order_relaxed);
                threadToSelect.store(omp_get_thread_num(), std::memory_order_relaxed);

                decreasePerThread[omp_get_thread_num()] = std::make_pair(v, farnessDecrease);
            }
        }

        const int thread = threadToSelect.load(std::memory_order_relaxed);

        if (thread != -1) { // There a convenient swap
            node v;
            WeightType farnessDecrease;
            std::tie(v, farnessDecrease) = decreasePerThread[thread];

            // Removing u from group
            groupFarness -= (farnessDecrease - increaseU);

            // Put v into the group
            group.insert(v);

            storeStatusFromThread(thread);
#ifdef NETWORKIT_SANITY_CHECKS
            runSanityCheck();
#endif // NETWORKIT_SANITY_CHECKS
            return true;
        }

        // Swaps with u did not improve the closeness enough, restore u into the group
        group.insert(u);
        const WeightType decreaseU = computeFarnessDecrease(u);
        tlx::unused(decreaseU);

        if (G->isWeighted())
            assert(std::abs(static_cast<double>(decreaseU - increaseU)) <= 1e-9);
        else
            assert(decreaseU == increaseU);

    } while (!candidatesToRemove.empty());

    return false;
}

template <class WeightType>
void GroupClosenessLocalSearchImpl<WeightType>::run() {

    if (runGrowShrink) {
        GroupClosenessGrowShrink gcgs(*G, group.begin(), group.end());
        gcgs.run();
        const auto newGroup = gcgs.groupMaxCloseness();
        group = std::unordered_set<node>{newGroup.begin(), newGroup.end()};
    }

    init();

    // Compute 1st and 2nd shortest distances and the vertices in the group that realize them
    computeDistances();

#ifdef NETWORKIT_SANITY_CHECKS
    runSanityCheck();
#endif // NETWORKIT_SANITY_CHECKS

    for (nIterations = 0; nIterations < maxIterations && findAndSwap(); ++nIterations) {
    }

    hasRun = true;
}

#ifdef NETWORKIT_SANITY_CHECKS
template <>
void GroupClosenessLocalSearchImpl<count>::runSanityCheck(node v, bool parallel) const {
    auto tmpGroup = group;
    if (v != none)
        tmpGroup.insert(v);

    const count n = G->upperNodeIdBound();
    std::vector<count> checkDistance(n, std::numeric_limits<count>::max());
    std::vector<count> checkDistance2(n, std::numeric_limits<count>::max());
    for (const node source : tmpGroup) {
        Traversal::BFSfrom(*G, source, [&](node u, count d) {
            if (d < checkDistance[u]) {
                assert(checkDistance2[u] >= checkDistance[u]);
                checkDistance2[u] = checkDistance[u];
                checkDistance[u] = d;
            } else if (d < checkDistance2[u])
                checkDistance2[u] = d;
        });
    }

    const auto &curDistance = parallel ? distanceGlobal[omp_get_thread_num()] : distance;
    const auto &curDistance2 = parallel ? distance2Global[omp_get_thread_num()] : distance2;

    G->forNodes([&](node u) {
        if (group.find(u) != group.end())
            return;
        assert(curDistance2[u] == checkDistance2[u]);
        assert(curDistance[u] == checkDistance[u]);
    });
}

template <>
void GroupClosenessLocalSearchImpl<edgeweight>::runSanityCheck(node v, bool parallel) const {
    auto tmpGroup = group;
    if (v != none)
        tmpGroup.insert(v);

    const count n = G->upperNodeIdBound();
    std::vector<edgeweight> checkDistance(n, std::numeric_limits<edgeweight>::max());
    std::vector<edgeweight> checkDistance2(n, std::numeric_limits<edgeweight>::max());
    for (const auto source : tmpGroup) {
        Traversal::DijkstraFrom(*G, source, [&](node u, edgeweight d) {
            if (d < checkDistance[u]) {
                checkDistance2[u] = checkDistance[u];
                checkDistance[u] = d;
            } else if (d < checkDistance2[u])
                checkDistance2[u] = d;
        });
    }

    const auto &curDistance = parallel ? distanceGlobal[omp_get_thread_num()] : distance;
    const auto &curDistance2 = parallel ? distance2Global[omp_get_thread_num()] : distance2;
    G->forNodes([&](node u) {
        if (group.find(u) != group.end())
            return;
        assert(std::abs(curDistance2[u] - checkDistance2[u]) <= 1e-9);
        assert(std::abs(curDistance[u] - checkDistance[u]) <= 1e-9);
    });
}
#endif // NETWORKIT_SANITY_CHECKS

} // namespace

GroupClosenessLocalSearch::GroupClosenessLocalSearch(const Graph &G, const std::vector<node> &group,
                                                     bool runGrowShrink, count maxIterations)
    : weighted(G.isWeighted()) {

    if (weighted)
        impl = std::make_unique<GroupClosenessLocalSearchImpl<edgeweight>>(
            G, group.begin(), group.end(), runGrowShrink, maxIterations);
    else
        impl = std::make_unique<GroupClosenessLocalSearchImpl<count>>(G, group.begin(), group.end(),
                                                                      runGrowShrink, maxIterations);
}

void GroupClosenessLocalSearch::run() {
    impl->run();
    hasRun = true;
}

std::vector<node> GroupClosenessLocalSearch::groupMaxCloseness() const {
    assureFinished();
    return {impl->group.begin(), impl->group.end()};
}

count GroupClosenessLocalSearch::numberOfIterations() const {
    assureFinished();
    return impl->nIterations;
}

} // namespace NetworKit
