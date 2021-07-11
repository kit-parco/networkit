/*
 * BarabasiAlbertGenerator.cpp
 *
 *  Created on: May 28, 2013
 *      Author: forigem, Manuel Penschuck <networkit@manuel.jetzt>
 */

#include <random>
#include <unordered_set>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/generators/BarabasiAlbertGenerator.hpp>

namespace NetworKit {
BarabasiAlbertGenerator::BarabasiAlbertGenerator(count k, count nMax, count n0, bool batagelj)
    : initGraph(0), k(k), nMax(nMax), batagelj(batagelj) {
    if (k > nMax)
        throw std::runtime_error("k (number of attachments per node) may not be larger than the "
                                 "number of nodes in the target graph (nMax)");
    if (n0 > nMax)
        throw std::runtime_error("n0 (number of initially connected nodes) may not be larger than "
                                 "the number of nodes in the target graph (nMax)");

    if (n0 < k) {
        if (n0 > 0) {
            WARN("given n0 is smaller than k, setting n0 = k");
        }
        this->n0 = k;
    } else {
        this->n0 = n0;
    }
}

BarabasiAlbertGenerator::BarabasiAlbertGenerator(count k, count nMax, const Graph &initGraph,
                                                 bool batagelj)
    : initGraph(initGraph), k(k), nMax(nMax), n0(0), batagelj(batagelj) {
    if (initGraph.numberOfNodes() != initGraph.upperNodeIdBound())
        throw std::runtime_error("initGraph is expected to have consecutive node ids");
    if (k > nMax)
        throw std::runtime_error("k (number of attachments per node) may not be larger than the "
                                 "number of nodes in the target graph (nMax)");
    if (initGraph.numberOfNodes() > nMax)
        throw std::runtime_error(
            "initialization graph cannot have more nodes than the target graph (nMax)");
    if (!batagelj && initGraph.numberOfNodes() < k) {
        throw std::runtime_error(
            "initialization graph for the original method needs at least k nodes");
    }
}

Graph BarabasiAlbertGenerator::generate() {
    if (!nMax)
        return Graph();

    if (batagelj) {
        return generateBatagelj();
    } else {
        return generateOriginal();
    }
}

Graph BarabasiAlbertGenerator::generateOriginal() {
    Graph G(nMax);
    if (n0 != 0) {
        // initialize the graph with n0 connected nodes
        for (count i = 1; i < n0; i++) {
            G.addEdge(i - 1, i);
        }
    } else {
        // initialize the graph with the edges from initGraph
        // and set n0 accordingly
        initGraph.forEdges([&G](node u, node v) { G.addEdge(u, v); });
        n0 = initGraph.upperNodeIdBound();
    }
    assert(G.numberOfNodes() >= k);

    Aux::SignalHandler handler;
    auto &gen = Aux::Random::getURNG();
    for (node u = n0; u < static_cast<node>(nMax); ++u) {
        std::uniform_int_distribution<uint64_t> indexDist{0, 2 * G.numberOfEdges()};

        std::unordered_set<node> targets;
        targets.reserve(k + 1);
        targets.insert(u);

        while (targets.size() - 1 < k) {
            auto randomIndex = indexDist(gen);

            for (node v : G.nodeRange()) {
                if (randomIndex <= G.degree(v)) { // found a node to connect to
                    targets.insert(v);
                    break;
                }
                randomIndex -= G.degree(v);
            }

            handler.assureRunning();
        }

        targets.erase(u);

        G.preallocateUndirected(u, k);
        for (node x : targets) {
            G.addEdge(u, x);
        }
    }

    G.shrinkToFit();
    return G;
}

Graph BarabasiAlbertGenerator::generateBatagelj() {
    const node n = nMax;

    // Temporarily stored edges in edge list M to allow fast  random access.
    // Degrees are only computed to accelerate graph building later on.
    // TODO: Once we've a fast GraphBuilder remove degree and migrate to GraphBuilder
    std::vector<node> M;
    std::vector<count> degree(n, 0);

    auto addEdge = [&](node u, node v) {
        M.push_back(u);
        M.push_back(v);
        degree[u]++;
        degree[v]++;
    };

    // copy seed graph into M
    if (initGraph.numberOfNodes() == 0) {
        M.reserve(2 * n0 + 2 * (n - n0) * k);

        // initialize n0 connected nodes
        for (index v = 0; v < n0 - 1; ++v) {
            addEdge(v, v + 1);
        }
        addEdge(0, n0 - 1);
    } else {
        M.reserve(2 * initGraph.numberOfEdges() + 2 * (n - initGraph.numberOfNodes()) * k);

        initGraph.forEdges([&](node u, node v) { addEdge(u, v); });
        n0 = initGraph.numberOfNodes();
    }

    // for each of the remaining nodes [n0, n), we draw k random DIFFERENT neighbors
    auto &gen = Aux::Random::getURNG();
    Aux::SignalHandler handler;
    for (index v = n0; v < n; ++v) {
        // If we were to update the range in the next loop, the additionally available nodes
        // would only lead to self-loops are multi-edges.
        std::uniform_int_distribution<size_t> distr{0, M.size() - 1};
        auto firstNeighbor = M.size() + 1;

        for (index i = 0; i < k; ++i) {
            // let's sample a new neighbor and repeat if we're already connected to it
            while (true) {
                const auto randomIndex = distr(gen);
                const auto newNeighbor = M[randomIndex];

                // the last 2*(i-1) positions contain all edges incident to v in to format
                //  Even  Odd   Even  Odd   Even  Odd
                // | v | Neigh | v | Neigh | v | Neigh ...
                // Hence, we need to compare the new neighbor to the previous (i-1) odd positions
                bool alreadyIncident = false;
                for (auto j = firstNeighbor; j < M.size(); j += 2) {
                    assert(M[j] != v); // ensure that we're not off by 1

                    if (M[j] == newNeighbor) {
                        alreadyIncident = true;
                        break;
                    }
                }

                if (!alreadyIncident) {
                    addEdge(v, newNeighbor);
                    break;
                }

                handler.assureRunning();
            }
        }
    }

    Graph G(nMax);

    for (node u = 0; u < n; ++u)
        G.preallocateUndirected(u, degree[u]);

    for (size_t i = 0; i < M.size(); i += 2)
        G.addEdge(M[i], M[i + 1]);

    return G;
}

} /* namespace NetworKit */
