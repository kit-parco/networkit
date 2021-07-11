// no-networkit-format
/*
 * GroupCloseness.cpp
 *
 *  Created on: 03.10.2016
 *      Author: elisabetta bergamini
 */

#include <networkit/centrality/GroupCloseness.hpp>
#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/centrality/TopCloseness.hpp>

#include <atomic>
#include <memory>
#include <omp.h>
#include <queue>

namespace NetworKit {

GroupCloseness::GroupCloseness(const Graph &G, count k, count H)
    : G(G), k(k), H(H) {}

edgeweight GroupCloseness::computeImprovement(node u, count n, Graph &G,
                                              count h) {
    // computes the marginal gain due to adding u to S
    std::vector<count> d1(n);
    G.forNodes([&](node v) { d1[v] = d[v]; });

    d1[u] = 0;
    count improvement = d[u]; // old distance of u
    std::queue<node> Q;
    Q.push(u);

    index level = 0;
    while (!Q.empty() && (h == 0 || level <= h)) {
        node v = Q.front();
        Q.pop();
        level = d1[v];
        G.forNeighborsOf(v, [&](node w) {
            if (d1[w] > d1[v] + 1) {
                d1[w] = d1[v] + 1;
                improvement += d[w] - d1[w];
                Q.push(w);
            }
        });
    }
    return improvement;
}

std::vector<count> GroupCloseness::newDistances(node u, count n, Graph &G,
                                                count h) {
    std::vector<count> d1(n);
    G.forNodes([&](node v) { d1[v] = d[v]; });

    d1[u] = 0;
    count improvement = d[u]; // old distance of u
    std::queue<node> Q;
    Q.push(u);

    index level = 0;
    while (!Q.empty() && (h == 0 || level <= h)) {
        node v = Q.front();
        Q.pop();
        level = d1[v];
        G.forNeighborsOf(v, [&](node w) {
            if (d1[w] > d1[v] + 1) {
                d1[w] = d1[v] + 1;
                improvement += d[w] - d1[w];
                Q.push(w);
            }
        });
    }
    return d1;
}

bool pairCompare(const std::pair<node, count> &firstElem,
                 const std::pair<node, count> &secondElem) {
    return firstElem.second > secondElem.second;
}

void GroupCloseness::run() {
    count n = G.upperNodeIdBound();
    node top = 0;
    iters = 0;
    std::vector<bool> visited(n, false);
    std::vector<node> pred(n);
    std::vector<count> distances(n);
    D.clear();
    D.resize(n, 0);
    std::vector<std::pair<node, count>> degPerNode(n);

    // compute degrees per node
    G.parallelForNodes([&](node v) {
        D[v] = G.degree(v);
        degPerNode[v] = std::make_pair(v, D[v]);
    });
    omp_lock_t lock;
    omp_init_lock(&lock);
    // sort by degree (in descending order) and retrieve max and argmax
    std::sort(degPerNode.begin(), degPerNode.end(), pairCompare);
    node nodeMaxDeg = degPerNode[0].first;

    if (H == 0) {
        TopCloseness topcc(G, 1, true, false);
        topcc.run();
        top = topcc.topkNodesList()[0];
    } else {
        top = nodeMaxDeg;
    }

    // first, we store the distances between each node and the top node
    d.clear();
    d.resize(n);
    BFS bfs(G, top);
    bfs.run();

    G.parallelForNodes([&](node v) { d[v] = bfs.distance(v); });

    // get max distance
    maxD = 0;
    count sumD = 0;
    // TODO: actually, we could have more generic parallel reduction iterators in
    // the Graph class
    G.forNodes([&](node v) {
        if (d[v] > maxD) {
            maxD = d[v];
        }
        sumD += d[v];
    });
    INFO("maxD = ", maxD);

    // init S
    S.clear();
    S.resize(k, 0);
    S[0] = top;

    std::vector<int64_t> prevBound(n, 0);
    d1.clear();
    count currentImpr = sumD + 1; // TODO change
    count maxNode = 0;

    std::vector<count> S2(n, sumD);
    S2[top] = 0;
    std::vector<int64_t> prios(n);

    // loop to find k group members
    for (index i = 1; i < k; i++) {
        DEBUG("k = ", i);
        G.parallelForNodes([&](node v) { prios[v] = -prevBound[v]; });
        // Aux::BucketPQ Q(prios, currentImpr + 1);
        Aux::BucketPQ Q(n, -currentImpr - 1, 0);
        G.forNodes([&](node v) { Q.insert(prios[v], v); });
        currentImpr = 0;
        maxNode = 0;
        d1.resize(G.upperNodeIdBound());

        std::atomic<bool> toInterrupt{false};
#pragma omp parallel // Shared variables:
        // cc: synchronized write, read leads to a positive race condition;
        // Q: fully synchronized;
        {
            while (!toInterrupt.load(std::memory_order_relaxed)) {
                omp_set_lock(&lock);
                if (Q.empty()) {
                    omp_unset_lock(&lock);
                    toInterrupt.store(true, std::memory_order_relaxed);
                    break;
                }

                auto topPair = Q.extractMin();
                node v = topPair.second;

                omp_unset_lock(&lock);
                INFO("Extracted node ", v, " with prio ", prevBound[v]);
                if (i > 1 && prevBound[v] <= static_cast<int64_t>(currentImpr)) {
                    INFO("Interrupting! currentImpr = ", currentImpr,
                         ", previous bound = ", prevBound[v]);
                    toInterrupt.store(true, std::memory_order_relaxed);
                    break;
                }
                if (D[v] > 1 && !(d[v] == 1 && D[v] == 2) && d[v] > 0 &&
                    (i == 1 || prevBound[v] > static_cast<int64_t>(currentImpr))) {
                    count imp = computeImprovement(v, n, G, H);
                    omp_set_lock(&lock);
                    if (imp > currentImpr) {
                        currentImpr = imp;
                        maxNode = v;
                    }
                    omp_unset_lock(&lock);
                    INFO("New bound for v = ", imp);
                    prevBound[v] = imp; // TODO use later
                } else {
                    prevBound[v] = 0;
                }
            }
        }
        S[i] = maxNode;

        d1 = newDistances(S[i], n, G, 0);
        G.parallelForNodes([&](node v) { d[v] = d1[v]; });
    }

    hasRun = true;
}

double GroupCloseness::computeFarness(std::vector<node> S, count H) {
    // we run a BFS from S up to distance H (if H > 0) and sum the distances
    double farness = 0;
    count k = S.size();
    std::vector<double> d1(G.upperNodeIdBound(), 0);
    std::vector<bool> visited(G.upperNodeIdBound(), false);
    std::queue<node> Q;

    for (node i = 0; i < k; i++) {
        Q.push(S[i]);
        visited[S[i]] = true;
    }

    while (!Q.empty()) {
        node u = Q.front();
        Q.pop();
        if (H > 0 && d1[u] > H) {
            break;
        }
        farness += d1[u];
        G.forNeighborsOf(u, [&](node w) {
            if (!visited[w]) {
                visited[w] = true;
                d1[w] = d1[u] + 1;
                Q.push(w);
            }
        });
    }
    return farness;
}

} /* namespace NetworKit */
