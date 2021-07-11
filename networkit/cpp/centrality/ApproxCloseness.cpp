// no-networkit-format
/*
 * ApproxCloseness.cpp
 *
 *  Created on: Dec 8, 2015
 *      Author: Sarah Lutteropp (uwcwa@student.kit.edu) and Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <networkit/centrality/ApproxCloseness.hpp>
#include <networkit/auxiliary/PrioQueue.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <cassert>
#include <queue>

namespace NetworKit {

struct ListEntry {
    node node_val;
    edgeweight dist_val;
};

ApproxCloseness::ApproxCloseness(const Graph& G, count nSamples, double epsilon, bool normalized, CLOSENESS_TYPE type) : Centrality(G, normalized), nSamples(nSamples), epsilon(epsilon), type(type) {
    assert(nSamples > 0 && epsilon >= 0);
}

void ApproxCloseness::run() {
    if (nSamples > G.numberOfNodes()) {
        WARN("Number of samples higher than the number of nodes. Setting number of samples to number of nodes");
        nSamples = G.numberOfNodes();
    }
    if (G.isDirected()) {
        switch (type) {
            case OUTBOUND:
                estimateClosenessForDirectedGraph(true);
                break;
            case INBOUND:
                estimateClosenessForDirectedGraph(false);
                break;
            case SUM:
            {
                estimateClosenessForDirectedGraph(true);
                std::vector<double> outbound = scoreData;
                estimateClosenessForDirectedGraph(false);
                G.parallelForNodes([&](node u) {
                    scoreData[u] += outbound[u];
                });
                break;
            }
            default:
                break;
        }

        G.parallelForNodes([&](node u) {
            if (fabs(scoreData[u]) > 1e-9) {
                scoreData[u] = 1/scoreData[u];
            }

        });

    } else {
        estimateClosenessForUndirectedGraph();
        G.parallelForNodes([&](node u) {
            scoreData[u] = normalized? (static_cast<double>(G.numberOfNodes()-1)) / scoreData[u] : 1 / scoreData[u];
        });
    }

    hasRun = true;
}

void ApproxCloseness::estimateClosenessForUndirectedGraph() {
    std::vector<node> sampledNodes = GraphTools::randomNodes(G, nSamples);

    LCSum = std::vector<double>(G.upperNodeIdBound());
    LCNum = std::vector<count>(G.upperNodeIdBound());
    LCSumSQ = std::vector<double>(G.upperNodeIdBound());
    HCSum = std::vector<double>(G.upperNodeIdBound());
    HCSumSQErr = std::vector<double>(G.upperNodeIdBound());
    HSum = std::vector<double>(G.upperNodeIdBound());
    HNum = std::vector<count>(G.upperNodeIdBound());
    scoreData = std::vector<double>(G.upperNodeIdBound());
    SQErrEst = std::vector<double> (G.upperNodeIdBound());

    std::vector<node> pivot(G.upperNodeIdBound());
    std::vector<edgeweight> delta(G.upperNodeIdBound());
    computeClosestPivot(sampledNodes, pivot, delta);

    for (count i = 0; i < nSamples; ++i) {
        runOnPivot(i, pivot, delta, sampledNodes);
    }

    G.parallelForNodes([&](node u) {
        if (sampledNodes[pivot[u]] != u) { // exclude sampled nodes
            const double LNum = static_cast<double>(G.numberOfNodes() - 1 - HNum[u] - nSamples + LCNum[u]);
            count HCNum = nSamples - LCNum[u];

            bool includeHCTerm = true;
            if (HCNum == 0) includeHCTerm = false;

            double p = static_cast<double>(LCNum[u]) / LNum;
            scoreData[u] = HSum[u] + HCSum[u] + LCSum[u] / p;
            double LCSqAvg = (LCSum[u] / static_cast<double>(LCNum[u]))
                             * (LCSum[u] / static_cast<double>(LCNum[u]));
            if (includeHCTerm) {
                SQErrEst[u] = 1.0 / static_cast<double>(LCNum[u])
                                  * (LCSumSQ[u] / static_cast<double>(LCNum[u]) - LCSqAvg) * LNum
                              + HCSumSQErr[u] / static_cast<double>(HCNum * HNum[u]);
            } else {
                SQErrEst[u] = 1.0 / static_cast<double>(LCNum[u])
                              * (LCSumSQ[u] / static_cast<double>(LCNum[u]) - LCSqAvg) * LNum;
            }

        } else {
            SQErrEst[u] = 0.0;
        }
    });
}

void ApproxCloseness::estimateClosenessForDirectedGraph(bool outbound) {
    if (G.isWeighted()) {
        computeClosenessForDirectedWeightedGraph(outbound);
    } else {
        computeClosenessForDirectedUnweightedGraph(outbound);
    }
}

void ApproxCloseness::computeClosenessForDirectedWeightedGraph(bool outbound) {
    R = std::vector<double>(G.upperNodeIdBound());
    scoreData = std::vector<double>(G.upperNodeIdBound());

    unsigned int t = 0;
    std::vector<bool> mark(G.upperNodeIdBound(), false);
    std::vector<unsigned int> count(G.upperNodeIdBound(), 0);
    std::vector<long> T(G.upperNodeIdBound(), 0);
    std::vector<edgeweight> distSum(G.upperNodeIdBound(), 0);
    std::vector<edgeweight> dist(G.upperNodeIdBound(), infDist);
    std::vector<index> round(G.upperNodeIdBound(), 0);
    Aux::PrioQueue<edgeweight, node> pq(dist.size());

    G.forNodesInRandomOrder([&](node u) {
        t++;
        mark[u] = true;
        // Perform pruned Dijkstra from u on G^T (or on G, if inbound)
        dist[u] = 0.0;
        round[u] = t;
        pq.insert(dist[u], u);
        while (!pq.empty()) {
            node v = pq.extractMin().second;
            if (dist[v] == infDist) break;
            if (count[v] < nSamples) {  // only continue if count[v] < k
                if (u != v) {
                    distSum[v] += dist[v];
                    count[v]++;
                    if (count[v] == nSamples) {
                        T[v] = t;
                        if (mark[v]) T[v] = t - 1;
                    }
                }

                // Dijkstra part
                if (outbound) {
                    G.forInNeighborsOf(v, [&](node v2, edgeweight w) {
                        if (round[v2] < t || dist[v] + w < dist[v2]) {
                            dist[v2] = dist[v] + w;
                            pq.changeKey(dist[v2], v2);
                            round[v2] = t;
                        }
                    });
                } else {
                    G.forNeighborsOf(v, [&](node v2, edgeweight w) {
                        if (round[v2] < t || dist[v] + w < dist[v2]) {
                            dist[v2] = dist[v] + w;
                            pq.changeKey(dist[v2], v2);
                            round[v2] = t;
                        }
                    });
                }
            }
        }
    });

    G.parallelForNodes([&](node v) {
        if (count[v] == 0) {
            scoreData[v] = 0;
        } else {
            scoreData[v] = distSum[v] / (double) count[v];
        }
        if (count[v] < nSamples) {
            R[v] = static_cast<double>(count[v]);
        } else {
            R[v] = 1 + static_cast<double>((nSamples - 1) * (G.numberOfNodes() - 2)) / static_cast<double>(T[v] - 1);
        }
    });
}

void ApproxCloseness::computeClosenessForDirectedUnweightedGraph(bool outbound) {
    R = std::vector<double>(G.upperNodeIdBound());
    scoreData = std::vector<double>(G.upperNodeIdBound());

    unsigned int t = 0;
    std::vector<bool> mark(G.upperNodeIdBound(), false);
    std::vector<unsigned int> count(G.upperNodeIdBound(), 0);
    std::vector<long> T(G.upperNodeIdBound(), 0);
    std::vector<edgeweight> distSum(G.upperNodeIdBound(), 0);
    std::vector<edgeweight> dist(G.upperNodeIdBound(), infDist);
    std::vector<index> round(G.upperNodeIdBound(), 0);

    G.forNodesInRandomOrder([&](node u) {
        t++;
        mark[u] = true;
        // Perform pruned BFS from u on G^T (or on G, if inbound)
        dist[u] = 0.0;
        round[u] = t;

        std::queue<node> q;
        q.push(u);
        while (!q.empty()) {
            node v = q.front(); q.pop();
            if (count[v] < nSamples) {  // only continue if count[v] < k
                if (u != v) {
                    distSum[v] += dist[v];
                    count[v]++;
                    if (count[v] == nSamples) {
                        T[v] = t;
                        if (mark[v]) T[v] = t - 1;
                    }
                }

                // BFS part
                if (outbound) {
                    G.forInNeighborsOf(v, [&](node v2, edgeweight) {
                        if (round[v2] < t) {
                            dist[v2] = dist[v] + 1;
                            q.push(v2);
                            round[v2] = t;
                        }
                    });
                } else {
                    G.forNeighborsOf(v, [&](node v2, edgeweight) {
                        if (round[v2] < t) {
                            dist[v2] = dist[v] + 1;
                            q.push(v2);
                            round[v2] = t;
                        }
                    });
                }
            }
        }

    });

    G.parallelForNodes([&](node v) {
        if (count[v] == 0) {
            scoreData[v] = 0;
        } else {
            scoreData[v] = distSum[v] / (double) count[v];
        }
        if (count[v] < nSamples) {
            R[v] = count[v];
        } else {
            R[v] = 1 + static_cast<double>((nSamples - 1) * (G.numberOfNodes() - 2)) / static_cast<double>(T[v] - 1);
        }
    });
}

void ApproxCloseness::computeClosestPivot(const std::vector<node> &samples, std::vector<node> &pivot, std::vector<edgeweight> &delta) {
    std::fill(delta.begin(), delta.end(), infDist);

    Aux::PrioQueue<edgeweight, node> pq(delta.size());
    for (index i = 0; i < samples.size(); ++i) {
        delta[samples[i]] = 0.0; // distance to closest pivot is 0 for pivot itself
        pivot[samples[i]] = i; // sample node is its own pivot
        pq.insert(0.0, samples[i]);
    }
    while (!pq.empty()) {
        node u = pq.extractMin().second;
        G.forNeighborsOf(u, [&](node v, edgeweight w) {
            if (delta[u] + w < delta[v]) {
                delta[v] = delta[u] + w;
                pivot[v] = pivot[u];
                pq.changeKey(delta[v], v);
            }
        });
    }
}

void ApproxCloseness::runOnPivot(index i, const std::vector<node> &pivot, const std::vector<edgeweight> &delta, const std::vector<node> &samples) {
    std::vector<edgeweight> pivotDist;
    std::vector<node> order;
    orderNodesByIncreasingDistance(samples[i], order, pivotDist);

    std::vector<node> last(samples.size(), G.upperNodeIdBound());
    std::vector<edgeweight> dist(samples.size());

    std::vector<std::vector<ListEntry>> list(samples.size());
    std::vector<std::vector<node>> nodes(G.upperNodeIdBound());

    std::vector<double> thresh(order.size());
    std::vector<double> bin(order.size());
    std::vector<count> count_vec(order.size());
    thresh[0] = 0;
    size_t curt = 0;
    size_t t = 0;

    for (node u : order) {
        if (pivotDist[u] == infDist) break; // all remaining nodes in the queue have infinite distance -> we are done.
        edgeweight d = pivotDist[u];
        scoreData[samples[i]] += d;
        if (samples[pivot[u]] == u) { // u belongs to the sampled nodes, compute its score exactly
            index j = pivot[u];
            last[j] = i;
            dist[j] = d;
            for (const ListEntry &z : list[j]) {
                if (epsilon != 0.0 && d > delta[z.node_val] / epsilon) {
                    HCSum[z.node_val] += z.dist_val;
                    HCSumSQErr[z.node_val] += (z.dist_val - d) * (z.dist_val - d);
                } else {
                    LCSum[z.node_val] += z.dist_val;
                    LCNum[z.node_val]++;
                    LCSumSQ[z.node_val] += z.dist_val * z.dist_val;
                }
            }
            list[j].clear();
        } else { // u is another node, estimate its score
            if (epsilon == 0 || (d <= delta[u] * (1.0 / epsilon - 1.0)) || ((last[pivot[u]] == i) && (dist[pivot[u]] <= delta[u] / epsilon))) {
                LCSum[u] += d;
                LCNum[u]++;
                LCSumSQ[u] += d * d;
            } else {
                list[pivot[u]].push_back({u, pivotDist[u]});
            }

            if (pivot[u] == i) {
                if (epsilon == 0.0 || fabs(thresh[t] - d / epsilon) < 1e-9) {
                    nodes[t].push_back(u);
                } else {
                    t++;
                    thresh[t] = (d / epsilon);
                    nodes[t].clear();
                    nodes[t].push_back(u);
                    bin[t] = 0;
                    count_vec[t] = 0;
                }
            }

            while (curt < t && d > thresh[curt + 1]) curt++;
            if (d > thresh[curt]) {
                bin[curt] += d;
                count_vec[curt]++;
            }
        }
    }

    // Compute tail sums for nodes for which c is pivot
    double tailsum = 0;
    count tailnum = 0;
    while (t > 0) {
        tailsum += bin[t];
        tailnum += count_vec[t];
        for (node u : nodes[t]) {
            HSum[u] = tailsum;
            HNum[u] = tailnum;
        }
        t--;
    }
}

void ApproxCloseness::orderNodesByIncreasingDistance(node c, std::vector<node> &order, std::vector<edgeweight> &pivotDist) {
    pivotDist = std::vector<edgeweight>(G.upperNodeIdBound(), infDist);
    pivotDist[c] = 0.0;
    order = std::vector<node>(G.numberOfNodes());

    if (G.isWeighted()) { // use Dijkstra
        Aux::PrioQueue<edgeweight, node> pq(pivotDist.size());
        pq.insert(0.0, c);
        index idx = 0;
        while (!pq.empty()) {
            node u = pq.extractMin().second;
            order[idx++] = u;
            G.forNeighborsOf(u, [&](node v, edgeweight w) {
                if (pivotDist[u] + w < pivotDist[v]) {
                    pivotDist[v] = pivotDist[u] + w;
                    pq.changeKey(pivotDist[v], v);
                }
            });
        }

        assert(idx == G.numberOfNodes());
    } else { // use BFS to compute distance from pivot and the respective order
        std::queue<node> q;
        q.push(c);
        pivotDist[c] = 0;
        std::vector<bool> visited(G.upperNodeIdBound(), false);
        visited[c] = true;
        index idx = 0;
        while (!q.empty()) {
            node u = q.front();
            q.pop();
            order[idx++] = u; // TODO: Is this correct?
            G.forNeighborsOf(u, [&](node v) {
                if (!visited[v]) {
                    pivotDist[v] = pivotDist[u] + 1;
                    q.push(v);
                    visited[v] = true;
                }
            });
        }
    }
}

double ApproxCloseness::maximum() {
    return 1. / static_cast<double>(G.numberOfNodes() - 1);
}

std::vector<double> ApproxCloseness::getSquareErrorEstimates() {
    assureFinished();
    return SQErrEst;
}

} /* namespace NetworKit */
