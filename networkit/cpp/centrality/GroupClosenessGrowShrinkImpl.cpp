/*
 *  GroupClosenessGrowShrinkImpl.cpp
 *
 *  Created on: 19.12.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#include <omp.h>
#include <utility>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>

#include "GroupClosenessGrowShrinkImpl.hpp"

namespace NetworKit {
namespace GroupClosenessGrowShrinkDetails {
template <class WeightType>
GroupClosenessGrowShrinkImpl<WeightType>::GroupClosenessGrowShrinkImpl(const Graph &G,
                                                                       std::vector<node> initGroup,
                                                                       bool extended,
                                                                       count insertions,
                                                                       count maxIterations)
    : G(&G), group(std::move(initGroup)), extended(extended),
      insertions(insertions == 0 ? computeConsecutiveInsertions(G, group.size()) : insertions),
      maxIterations(maxIterations), heap(CompareDistance(distance)),
      heap_(CompareDistance(distance_)) {

    if (G.isDirected())
        std::runtime_error("Error, this algorithm does not support directed graphs.");

#if __AVX2__
    INFO("AVX2 is available.");
#endif
}

template <class WeightType>
count GroupClosenessGrowShrinkImpl<WeightType>::computeConsecutiveInsertions(const Graph &G,
                                                                             size_t groupSize) {
    if (groupSize == 0)
        throw std::runtime_error("Error, empty group.");

    Diameter diam(G, DiameterAlgo::estimatedRange, 0.1);
    diam.run();
    const auto diamOverK =
        static_cast<double>(diam.getDiameter().second) / std::sqrt(static_cast<double>(groupSize));
    return std::max(1., .5 + diamOverK);
}

template <class WeightType>
void GroupClosenessGrowShrinkImpl<WeightType>::init() {
    const count n = G->upperNodeIdBound();
    distance.assign(n, 0);
    distance_.assign(n, 0);
    visited.assign(n, 0);
    sumOfMins.assign(n, 0);

    heap.reserve(n);
    if (G->isWeighted())
        heap_.reserve(n);

    nearest.assign(n, none);
    nearest_.assign(n, none);
    stack.assign(n, 0);

    intDistributions.resize(omp_get_max_threads());
    urngs.reserve(omp_get_max_threads());
    for (index i = 0; i < static_cast<index>(omp_get_max_threads()); ++i)
        urngs.emplace_back(Aux::Random::getURNG());

    for (size_t i = 0; i < group.size(); ++i)
        idxMap.emplace(group[i], i);

#ifdef __AVX2__
    randVec.resize(n);
#else
    randVec.resize(K * n);
#endif // __AVX2__

    totalSwaps = 0;
}

template <class WeightType>
void GroupClosenessGrowShrinkImpl<WeightType>::initRandomVec() {
    G->parallelForNodes([&](node u) {
        // Avoid to generate numbers for nodes in the group
        if (distance[u]) {
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

template <class WeightType>
void GroupClosenessGrowShrinkImpl<WeightType>::dijkstra() {
    auto &pq = G->isWeighted() ? heap_ : heap;
    assert(!pq.empty());

    do {
        const node x = pq.extract_top();

        G->forNeighborsOf(x, [&](node y, edgeweight w) {
            if (!visited[y] || distance_[y] > distance_[x] + static_cast<WeightType>(w)) {
                distance_[y] = distance_[x] + w;
                nearest_[y] = nearest_[x];
                pq.update(y);
                visited[y] = true;
            }
        });

    } while (!pq.empty());
}

template <class WeightType>
bool GroupClosenessGrowShrinkImpl<WeightType>::findAndSwap() {
    // Nodes that have been inserted/removed from the group
    std::vector<node> nodesInserted, nodesRemoved;
    nodesInserted.reserve(insertions);
    nodesRemoved.reserve(insertions);

    // Total farness decrement/increment after having inserted/removed the nodes from the group
    WeightType totalDecrement{0}, totalIncrement{0};

    for (count i = 0; i < insertions; ++i) {
        // Find node with highest farness decrement by estimating the size
        // of its BFS DAG
        const node v = estimateHighestDecrement();
        nodesInserted.push_back(v);

        // Put v in the group
        idxMap[v] = nextIdx.front();
        nextIdx.pop();

        // Update farness decrement
        totalDecrement += computeFarnessDecrement(v);
    }

    for (count i = 0; i < insertions; ++i) {
        // Update the farness decrement for each node in the group
        computeFarnessIncrement();

        // Compute node with lowest farness increment
        auto u = idxMap.begin()->first;
        auto minimumIncrement = increment[idxMap.begin()->second];

        std::for_each(++idxMap.begin(), idxMap.end(), [&](const std::pair<node, index> &entry) {
            if (increment[entry.second] < minimumIncrement) {
                u = entry.first;
                minimumIncrement = increment[entry.second];
            }
        });

        // Update farness increment
        totalIncrement += minimumIncrement;

        // Removing u from group
        nextIdx.push(idxMap.at(u));
        idxMap.erase(u);
        nodesRemoved.push_back(u);

        // Updating (1st) shortest distances and setting 2nd shortest distances to infinity
        G->forNodes([&](node x) {
            if (nearest[x] == u) {
                nearest[x] = nearest_[x];
                distance[x] = distance_[x];
                nearest_[x] = none;
                distance_[x] = infdistance;
            } else if (nearest_[x] == u) {
                nearest_[x] = none;
                distance_[x] = infdistance;
            }
        });

        auto &pq = G->isWeighted() ? heap_ : heap;

        // Update the distance of y
        const auto process = [&](node x, node y, WeightType w) -> void {
            if (nearest[x] != nearest[y]) {
                if (distance_[y] > distance[x] + w) {
                    distance_[y] = distance[x] + w;
                    nearest_[y] = nearest[x];
                    pq.update(y);
                }
                // Checking 2nd distanceance
            } else if (distance_[x] < infdistance && distance_[y] > distance_[x] + w) {
                distance_[y] = distance_[x] + w;
                nearest_[y] = nearest_[x];
                pq.update(y);
                assert(nearest_[x] != nearest[y]);
            }
        };

        G->forEdges([&](node x, node y, edgeweight w) {
            process(x, y, static_cast<WeightType>(w));
            process(y, x, static_cast<WeightType>(w));
        });

        do {
            const node x = pq.extract_top();
            G->forNeighborsOf(x, [&](node y, edgeweight w) {
                if (nearest[y] != nearest_[x]
                    && distance_[y] > distance_[x] + static_cast<WeightType>(w)) {
                    distance_[y] = distance_[x] + w;
                    nearest_[y] = nearest_[x];
                    pq.update(y);
                }
            });
        } while (!pq.empty());
    }

    if (totalDecrement <= totalIncrement) {
        // Farness could not be decreased, restore original state
        do {
            idxMap[nodesRemoved.back()] = 0;
            nodesRemoved.pop_back();
        } while (!nodesRemoved.empty());

        do {
            idxMap.erase(nodesInserted.back());
            nodesInserted.pop_back();
        } while (!nodesInserted.empty());

        assert(idxMap.size() == group.size());

        return false;
    }

    return true;
}

template <class WeightType>
void GroupClosenessGrowShrinkImpl<WeightType>::bfsFromGroup() {
    assert(heap.empty());
    std::fill(visited.begin(), visited.end(), false);

    // Only needed by algorithm for unweighted graphs
    std::queue<node> q;

    for (const auto &entry : idxMap) {
        if (G->isWeighted())
            heap.push(entry.first);
        else
            q.push(entry.first);

        distance[entry.first] = 0;
        visited[entry.first] = true;
        nearest[entry.first] = entry.first;
    }

    do {
        const auto x = G->isWeighted() ? heap.extract_top() : extractQueueTop(q);

        G->forNeighborsOf(x, [&](node y, edgeweight w) {
            if (!visited[y] || (G->isWeighted() && distance[y] > distance[x] + w)) {
                distance[y] = distance[x] + static_cast<WeightType>(w);
                nearest[y] = nearest[x];
                visited[y] = true;
                if (G->isWeighted())
                    heap.update(y);
                else
                    q.push(y);
            }
        });
    } while (!(G->isWeighted() ? heap.empty() : q.empty()));
}

template <class WeightType>
void GroupClosenessGrowShrinkImpl<WeightType>::computeFarnessIncrement() {
    std::fill(increment.begin(), increment.end(), 0);
    G->forNodes([&](node x) { increment[idxMap.at(nearest[x])] += distance_[x] - distance[x]; });
}

template <class WeightType>
WeightType GroupClosenessGrowShrinkImpl<WeightType>::computeFarnessDecrement(node v) {
    distance_[v] = distance[v];
    nearest_[v] = nearest[v];
    distance[v] = WeightType{0};
    nearest[v] = v;

    // Only needed by the algorithm for unweighted graphs
    std::queue<node> q;

    if (G->isWeighted()) {
        heap.push(v);
    } else {
        q.push(v);
        std::fill(visited.begin(), visited.end(), false);
        visited[v] = true;
    }

    WeightType decr = G->isWeighted() ? 0 : distance[v];

    auto processWeighted = [&](node x, node y, WeightType w) -> void {
        if (distance[y] > distance[x] + w) {
            // Nearest nodes could have been already updated
            if (nearest[y] != v) {
                nearest_[y] = nearest[y];
                nearest[y] = v;
                distance_[y] = distance[y];
            }
            distance[y] = distance[x] + w;
            heap.update(y);
        } else if (nearest[x] == v && nearest[y] != v && distance_[y] > distance[x] + w) {
            distance_[y] = distance[x] + w;
            nearest_[y] = v;
            heap.update(y);
        } else if (nearest_[x] == v && nearest[y] != v && distance_[y] > distance_[x] + w) {
            distance_[y] = distance_[x] + w;
            nearest_[y] = v;
            heap.update(y);
        }
    };

    auto processUnweighted = [&](node x, node y) -> void {
        if (visited[y])
            return;

        if (distance[y] > distance[x] + 1) {
            distance_[y] = distance[y];
            nearest_[y] = nearest[y];
            decr += distance[y] - (distance[x] + 1);
            distance[y] = distance[x] + 1;
            nearest[y] = v;
            q.push(y);
        } else if (nearest[x] == v && distance_[y] > distance[x] + 1) {
            distance_[y] = distance[x] + 1;
            nearest_[y] = v;
            q.push(y);
        } else if (nearest_[x] == v && distance_[y] > distance_[x] + 1) {
            distance_[y] = distance_[x] + 1;
            nearest_[y] = v;
            q.push(y);
        }

        visited[y] = true;
    };

    do {
        const node x = G->isWeighted() ? heap.extract_top() : extractQueueTop(q);

        if (G->isWeighted() && nearest[x] == v)
            decr += distance_[x] - distance[x];

        G->forNeighborsOf(x, [&](node y, edgeweight w) {
            if (G->isWeighted())
                processWeighted(x, y, static_cast<WeightType>(w));
            else
                processUnweighted(x, y);
        });
    } while (!(G->isWeighted() ? heap.empty() : q.empty()));

    return decr;
}

template <class WeightType>
node GroupClosenessGrowShrinkImpl<WeightType>::estimateHighestDecrement() {
    assert(heap.empty());
    std::fill(visited.begin(), visited.end(), false);

    // Needed only by the algorithm for unweighted graphs
    std::queue<node> q;
    stackSize = 0;
    for (const auto &idx : idxMap) {
        if (G->isWeighted())
            heap.push(idx.first);
        else
            q.push(idx.first);
        visited[idx.first] = true;
    }

    // Recompute the DAG
    do {
        const node x = G->isWeighted() ? heap.extract_top() : extractQueueTop(q);

        bool leaf = true;
        G->forNeighborsOf(x, [&](node y, edgeweight w) {
            if (!visited[y] || (G->isWeighted() && distance[y] > distance[x] + w)) {
                visited[y] = true;
                if (G->isWeighted()) {
                    heap.update(y);
                } else {
                    leaf = false;
                    q.push(y);
                }
            }
        });

        if (distance[x] && (G->isWeighted() || (!leaf || distance[x] == WeightType{1})))
            stack[stackSize++] = x;

    } while (!(G->isWeighted() ? heap.empty() : q.empty()));

    // Generating all the necessary random numbers
    initRandomVec();

    const auto isCandidate = [&](const node x) -> bool {
        return G->isWeighted() ? distance[x]
                               : (distance[x] == WeightType{1} || (extended && distance[x] > 1));
    };

    // Do 16 packed repetitions at once
    for (size_t i = 0; i < stackSize; ++i) {
        const node x = stack[stackSize - 1 - i];
#ifdef __AVX2__
        // 16 randomly generated integers;
        __m256i &x1 = randVec[x].vec;
        // Pulling leaves
        G->forNeighborsOf(x, [&](node y, edgeweight w) {
            if (distance[y] == distance[x] + w) {
                const __m256i &y1 = randVec[y].vec;
                x1 = _mm256_min_epu16(x1, y1);
            }
        });
        *(__m256i *)(&randVec[x].items) = x1;
#else
        G->forNeighborsOf(x, [&](node y, edgeweight w) {
            if (distance[y] == distance[x] + w)
                for (index i = 0; i < K; ++i)
                    randVec[K * x + i] = std::min(randVec[K * x + i], randVec[K * y + i]);
        });
#endif // __AVX2__

        if (isCandidate(x)) {
            sumOfMins[x] = 0;
            for (index j = 0; j < K; ++j) {
#ifdef __AVX2__
                sumOfMins[x] += randVec[x].items[j];
#else
                sumOfMins[x] = randVec[K * x + j];
#endif // __AVX2__
            }

            if (!sumOfMins[x])
                sumOfMins[x] = 1;
        }
    }

    node v = none;
    float bestEstimate = -1.f;
    G->forNodes([&](node x) {
        if (isCandidate(x) && sumOfMins[x]) {
            const auto estimate =
                static_cast<float>(distance[x])
                * (static_cast<float>(K) / (static_cast<float>(sumOfMins[x]) / maxInt16) - 1.f);
            if (estimate > bestEstimate) {
                v = x;
                bestEstimate = estimate;
            }
        }
    });

    assert(v != none);
    return v;
}

template <class WeightType>
node GroupClosenessGrowShrinkImpl<WeightType>::extractQueueTop(std::queue<node> &q) {
    const node u = q.front();
    q.pop();
    return u;
}

template <class WeightType>
void GroupClosenessGrowShrinkImpl<WeightType>::run() {
    init();

    increment.assign(group.size() + insertions, 0);
    for (count i = 0; i < insertions; ++i)
        nextIdx.push(group.size() + i);

    // Compute 1st shortest distances
    bfsFromGroup();

    // Compute 2nd shortest distances
    // Step 1: find 'bundary' nodes (sources of Dijkstra)
    auto &pq = (G->isWeighted() ? heap_ : heap);
    std::fill(visited.begin(), visited.end(), false);
    const auto process = [&](node x, node y, WeightType w) -> void {
        if (!visited[y] || distance_[y] > distance[x] + w) {
            distance_[y] = distance[x] + w;
            visited[y] = true;
            nearest_[y] = nearest[x];
            pq.update(y);
        }
    };

    G->forEdges([&](node x, node y, edgeweight w) {
        if (nearest[x] != nearest[y]) {
            process(x, y, static_cast<WeightType>(w));
            process(y, x, static_cast<WeightType>(w));
        }
    });

    // Step 2: Explore rest of the graph with Dijkstra
    dijkstra();

    // The main algorithm starts here
    while (findAndSwap() && totalSwaps++ < maxIterations) {
    }
}

template <class WeightType>
std::vector<node> GroupClosenessGrowShrinkImpl<WeightType>::groupMaxCloseness() const {
    std::vector<node> result;
    result.reserve(group.size());

    for (auto &entry : idxMap)
        result.push_back(entry.first);

    return result;
}

template <class WeightType>
count GroupClosenessGrowShrinkImpl<WeightType>::numberOfIterations() const {
    return totalSwaps;
}

template class GroupClosenessGrowShrinkImpl<count>;
template class GroupClosenessGrowShrinkImpl<edgeweight>;
} // namespace GroupClosenessGrowShrinkDetails
} // namespace NetworKit
