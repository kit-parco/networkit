
#include <algorithm>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

namespace GraphTools {

count computeMaxDegree(const Graph &G, bool inDegree = false) {
    count result = 0;
#ifndef NETWORKIT_OMP2
#pragma omp parallel for reduction(max : result)
    for (omp_index u = 0; u < static_cast<omp_index>(G.upperNodeIdBound()); ++u) {
        result = std::max(result, inDegree ? G.degreeIn(u) : G.degreeOut(u));
    }
#else
    G.forNodes([&](const node u) {
        result = std::max(result, inDegree ? G.degreeIn(u) : G.degreeOut(u));
    });
#endif
    return result;
}

edgeweight computeMaxWeightedDegree(const Graph &G, bool inDegree = false) {
    edgeweight result = 0;
#ifndef NETWORKIT_OMP2
#pragma omp parallel for reduction(max : result)
    for (omp_index u = 0; u < static_cast<omp_index>(G.upperNodeIdBound()); ++u) {
        result = std::max(result, inDegree ? G.weightedDegreeIn(u) : G.weightedDegree(u));
    }
#else
    G.forNodes([&](const node u) {
        result = std::max(result, inDegree ? G.weightedDegreeIn(u) : G.weightedDegree(u));
    });
#endif
    return result;
}

count maxDegree(const Graph &G) {
    return computeMaxDegree(G);
}

count maxInDegree(const Graph &G) {
    return computeMaxDegree(G, true);
}

edgeweight maxWeightedDegree(const Graph &G) {
    return computeMaxWeightedDegree(G);
}

edgeweight maxWeightedInDegree(const Graph &G) {
    return computeMaxWeightedDegree(G, true);
}

node randomNode(const Graph &G) {
    if (!G.numberOfNodes())
        return none;

    auto &gen = Aux::Random::getURNG();
    std::uniform_int_distribution<node> distr{0, G.upperNodeIdBound() - 1};
    node v;

    do {
        // When there are many deleted nodes, we might call Aux::Random::integer
        // many times, and it is very expensive.
        v = distr(gen);
    } while (!G.hasNode(v));

    return v;
}

std::vector<node> randomNodes(const Graph &G, count n) {
    assert(n <= G.numberOfNodes());
    std::vector<node> selectedNodes;
    std::vector<bool> alreadySelected(G.numberOfNodes(), false);

    if (n == G.numberOfNodes()) {
        selectedNodes.insert(selectedNodes.begin(), G.nodeRange().begin(), G.nodeRange().end());
    } else if (n > G.numberOfNodes() / 2) { // in order to minimize the calls to randomNode
                                            // we randomize the ones that aren't pivot
                                            // if the are more to be selected than not-selected
        std::fill(alreadySelected.begin(), alreadySelected.end(), true);

        for (count i = 0; i < G.numberOfNodes() - n; ++i) { // we have to sample distinct nodes
            node v = GraphTools::randomNode(G);
            while (!alreadySelected[v]) {
                v = GraphTools::randomNode(G);
            }
            alreadySelected[v] = false;
        }

        for (const auto sample : G.nodeRange()) {
            if (alreadySelected[sample]) {
                selectedNodes.push_back(sample);
                if (selectedNodes.size() == n)
                    break;
            }
        }
    } else {
        for (count i = 0; i < n; ++i) { // we have to selected distinct nodes
            node v = GraphTools::randomNode(G);
            while (alreadySelected[v]) {
                v = GraphTools::randomNode(G);
            }
            selectedNodes.push_back(v);
            alreadySelected[v] = true;
        }
    }
    return selectedNodes;
}

std::pair<node, node> randomEdge(const Graph &G, bool uniformDistribution) {
    if (!G.numberOfEdges()) {
        throw std::runtime_error("Error: the graph has no edges!");
    }

    if (uniformDistribution) {
        /*
         * The simple idea here is to interpret all neighborhoods next to each other, resulting
         * in a virtual vector of size m. Then we draw a random index and return the edge.
         * For undirected edges, the vector has size 2m; but the idea remains. There is one minor
         * complication for undirected edges with self-loops: each edge {u,v} with u != v is stored
         * twice (once in the neighborhood of u, once in v) but a loop (u, u) is only stored once.
         * To equalize the probabilities we reject edges {u,v} with u > v and try again. This leads
         * to less than two expected trails in and is only done for undirected graphs with
         * self-loops.
         */

        do {
            const auto upper =
                G.isDirected() ? G.numberOfEdges() : 2 * G.numberOfEdges() - G.numberOfSelfLoops();
            auto idx = Aux::Random::index(upper);

            node u, v;

            if (idx > upper / 2) {
                // assuming degrees are somewhat distributed uniformly, it's better to start with
                // larger nodes for large indices. In this case we have to mirror the index:
                idx = (upper - 1) - idx;

                for (u = G.upperNodeIdBound() - 1; idx >= G.degree(u); --u) {
                    idx -= G.degree(u);
                }

                v = G.getIthNeighbor(u, G.degree(u) - 1 - idx);

            } else {
                for (u = 0; idx >= G.degree(u); ++u) {
                    assert(u < G.upperNodeIdBound());
                    idx -= G.degree(u);
                }

                v = G.getIthNeighbor(u, idx);
            }

            if (G.numberOfSelfLoops() && !G.isDirected() && u > v)
                // reject (see above)
                continue;

            return {u, v};
        } while (true);
    }

    node u; // we will return edge (u, v)

    // fast way, but not a uniform random edge!
    do {
        u = GraphTools::randomNode(G);
    } while (!G.degree(u));

    const auto v = GraphTools::randomNeighbor(G, u);

    return {u, v};
}

std::vector<std::pair<node, node>> randomEdges(const Graph &G, count nr) {
    if (!nr)
        return {};

    if (!G.numberOfEdges()) {
        throw std::runtime_error(
            "Graph has no edges to sample from. Add edges to the graph first.");
    }

    std::vector<std::pair<node, node>> edges;

    auto &gen = Aux::Random::getURNG();
    std::vector<count> outDeg(G.upperNodeIdBound());
    G.forNodes([&outDeg, &G](const node u) { outDeg[u] = G.degree(u); });

    std::discrete_distribution<count> distribution(outDeg.begin(), outDeg.end());

    for (index i = 0; i < nr; i++) {
        node u, v; // we will pick edge (u, v)
        if (G.isDirected()) {
            u = distribution(gen);
            // should always be the case as  without
            // edges should have probability 0
            assert(G.degree(u));
            v = randomNeighbor(G, u);
        } else {
            // self-loops which appear only once in the outEdge arrays
            // easiest way it to ignore edges (u, v) with u > v
            do {
                u = distribution(gen);
                // should always be the case as  without
                // edges should have probability 0
                assert(G.degree(u));
                v = randomNeighbor(G, u);
            } while (u > v);
        }
        edges.emplace_back(u, v);
    }

    return edges;
}

node randomNeighbor(const Graph &G, node u) {
    if (!G.degree(u))
        return none;

    return G.getIthNeighbor(u, Aux::Random::integer(G.degree(u) - 1));
}

std::pair<node, node> size(const Graph &G) noexcept {
    return {G.numberOfNodes(), G.numberOfEdges()};
}

double density(const Graph &G) noexcept {
    if (G.numberOfNodes() <= 1)
        return 0;
    const auto n = static_cast<double>(G.numberOfNodes());
    const auto m =
        static_cast<double>((G.numberOfEdges() - G.numberOfSelfLoops()) * (G.isDirected() ? 1 : 2));
    return m / (n * (n - 1));
}

double volume(const Graph &G) {
    double volume = G.totalEdgeWeight();

    if (!G.isDirected()) {
        volume *= 2.0;
    }
    return volume;
}

Graph copyNodes(const Graph &G) {
    Graph C(G.upperNodeIdBound(), G.isWeighted(), G.isDirected());
    for (node u = 0; u < G.upperNodeIdBound(); ++u) {
        if (!G.hasNode(u)) {
            C.removeNode(u);
        }
    }
    return C;
}

Graph subgraphFromNodes(const Graph &G, const std::unordered_set<node> &nodes) {
    return subgraphFromNodes(G, nodes.begin(), nodes.end(), false);
}

Graph subgraphFromNodes(const Graph &G, const std::unordered_set<node> &nodes,
                        bool includeOutNeighbors, bool includeInNeighbors) {
    return subgraphAndNeighborsFromNodes(G, nodes, includeOutNeighbors, includeInNeighbors);
}

Graph subgraphAndNeighborsFromNodes(const Graph &G, const std::unordered_set<node> &nodes,
                                    bool includeOutNeighbors, bool includeInNeighbors) {
    const auto neighbors = [&] {
        std::unordered_set<node> neighbors;

        if (!includeOutNeighbors && !includeInNeighbors)
            return neighbors;

        for (node u : nodes) {
            if (includeOutNeighbors)
                for (const node v : G.neighborRange(u))
                    neighbors.insert(v);

            if (includeInNeighbors)
                for (const node v : G.inNeighborRange(u))
                    neighbors.insert(v);
        }

        return neighbors;
    }();

    /*
     * returns one of three different relevance scores:
     * 2: Is in the original nodes set
     * 1: Is a relevant neighbor (i.e., respective include*Neighbor was set)
     * 0: Neither of both
     */
    auto isRelevantNode = [&](const node u) {
        if (nodes.find(u) != nodes.end())
            return 2;
        if (!neighbors.empty() && neighbors.find(u) != neighbors.end())
            return 1;
        return 0;
    };

    Graph S(G.upperNodeIdBound(), G.isWeighted(), G.isDirected());
    // delete all nodes that are not in the node set
    S.forNodes([&](node u) {
        if (!isRelevantNode(u)) {
            S.removeNode(u);
        }
    });

    G.forEdges([&](node u, node v, edgeweight w) {
        // only include edges if at least one endpoint is in nodes (relevance 2),
        // and the other is either in nodes or in neighbors (relevance >= 1)
        if (isRelevantNode(u) + isRelevantNode(v) > 2) {
            S.addEdge(u, v, w);
        }
    });

    return S;
}

Graph toUndirected(const Graph &G) {
    if (!G.isDirected()) {
        WARN("The graph is already undirected");
    }

    return Graph(G, G.isWeighted(), false);
}

Graph toUnweighted(const Graph &G) {
    if (!G.isWeighted()) {
        WARN("The graph is already unweighted");
    }

    return Graph(G, false, G.isDirected());
}

Graph toWeighted(const Graph &G) {
    if (G.isWeighted()) {
        WARN("The graph is already weighted");
    }

    return Graph(G, true, G.isDirected());
}

Graph transpose(const Graph &G) {
    if (!G.isDirected()) {
        throw std::runtime_error("The transpose of an undirected graph is "
                                 "identical to the original graph.");
    }

    Graph GTranspose(G.upperNodeIdBound(), G.isWeighted(), true);

    // prepare edge id storage if input has indexed edges
    if (G.hasEdgeIds()) {
        GTranspose.indexEdges();
    }

#pragma omp parallel for
    for (omp_index u = 0; u < static_cast<omp_index>(G.upperNodeIdBound()); ++u) {
        if (G.hasNode(u)) {
            GTranspose.preallocateDirected(u, G.degreeIn(u), G.degreeOut(u));

            G.forInEdgesOf(u, [&](node, node v, edgeweight w, edgeid id) {
                GTranspose.addPartialOutEdge(unsafe, u, v, w, id);
            });

            G.forEdgesOf(u, [&](node, node v, edgeweight w, edgeid id) {
                GTranspose.addPartialInEdge(unsafe, u, v, w, id);
            });

        } else {
#pragma omp critical
            GTranspose.removeNode(u);
        }
    }

    GTranspose.setEdgeCount(unsafe, G.numberOfEdges());
    GTranspose.setNumberOfSelfLoops(unsafe, G.numberOfSelfLoops());
    GTranspose.setUpperEdgeIdBound(unsafe, G.upperEdgeIdBound());
    assert(GTranspose.checkConsistency());

    return GTranspose;
}

void append(Graph &G, const Graph &G1) {
    std::unordered_map<node, node> nodeMap;
    G1.forNodes([&](node u) {
        node u_ = G.addNode();
        nodeMap[u] = u_;
    });

    if (G.isWeighted()) {
        G1.forEdges([&](node u, node v, edgeweight ew) { G.addEdge(nodeMap[u], nodeMap[v], ew); });
    } else {
        G1.forEdges([&](node u, node v) { G.addEdge(nodeMap[u], nodeMap[v]); });
    }
}

void merge(Graph &G, const Graph &G1) {
    if (G1.upperNodeIdBound() > G.upperNodeIdBound()) {
        count prevBound = G.upperNodeIdBound();
        for (node i = prevBound; i < G1.upperNodeIdBound(); ++i) {
            G.addNode();
        }
        for (node i = prevBound; i < G1.upperNodeIdBound(); ++i) {
            if (!G1.hasNode(i)) {
                G.removeNode(i);
            }
        }
    }

    for (node i = 0; i < G.upperNodeIdBound(); ++i) {
        if (!G.hasNode(i) && G1.hasNode(i)) {
            G.restoreNode(i);
        }
    }

    G1.forEdges([&](node u, node v, edgeweight w) {
        // naive implementation takes $O(m \cdot d)$ for $m$ edges and max. degree
        // $d$ in this graph
        if (!G.hasEdge(u, v)) {
            G.addEdge(u, v, w);
        }
    });
}

Graph getCompactedGraph(const Graph &graph, const std::unordered_map<node, node> &nodeIdMap) {
    return getRemappedGraph(graph, nodeIdMap.size(), [&](node u) {
        const auto it = nodeIdMap.find(u);
        assert(it != nodeIdMap.cend());
        return it->second;
    });
}

std::unordered_map<node, node> getContinuousNodeIds(const Graph &graph) {
    std::unordered_map<node, node> nodeIdMap;
    count continuousId = 0;
    auto addToMap = [&nodeIdMap, &continuousId](node v) {
        nodeIdMap.insert(std::make_pair(v, continuousId++));
    };
    graph.forNodes(addToMap);
    return nodeIdMap;
}
std::unordered_map<node, node> getRandomContinuousNodeIds(const Graph &graph) {
    std::unordered_map<node, node> nodeIdMap;
    std::vector<node> nodes;
    nodes.reserve(graph.numberOfNodes());

    graph.forNodes([&](node u) { nodes.push_back(u); });

    std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());

    count continuousId = 0;
    for (node v : nodes) {
        nodeIdMap.insert(std::make_pair(v, continuousId++));
    };

    return nodeIdMap;
}

std::vector<node> invertContinuousNodeIds(const std::unordered_map<node, node> &nodeIdMap,
                                          const Graph &G) {
    assert(nodeIdMap.size() == G.numberOfNodes());
    std::vector<node> invertedIdMap(G.numberOfNodes() + 1);
    // store upper node id bound
    invertedIdMap[G.numberOfNodes()] = G.upperNodeIdBound();
    // inverted node mapping
    for (const auto &x : nodeIdMap) {
        invertedIdMap[x.second] = x.first;
    }
    return invertedIdMap;
}

Graph restoreGraph(const std::vector<node> &invertedIdMap, const Graph &G) {
    // with the inverted id map and the compacted graph, generate the original graph again
    Graph Goriginal(invertedIdMap.back(), G.isWeighted(), G.isDirected());
    index current = 0;
    Goriginal.forNodes([&](node u) {
        if (invertedIdMap[current] == u) {
            G.forNeighborsOf(current, [&](node v) { Goriginal.addEdge(u, invertedIdMap[v]); });
            ++current;
        } else {
            Goriginal.removeNode(u);
        }
    });
    return Goriginal;
}

void sortEdgesByWeight(Graph &G, bool decreasing) {
    if (decreasing)
        G.sortEdges([](auto e1, auto e2) {
            if (e1.weight == e2.weight)
                return e1.v < e2.v;
            return e1.weight > e2.weight;
        });
    else
        G.sortEdges([](auto e1, auto e2) {
            if (e1.weight == e2.weight)
                return e1.v < e2.v;
            return e1.weight < e2.weight;
        });
}

} // namespace GraphTools
} // namespace NetworKit
