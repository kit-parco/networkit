// no-networkit-format
/*
 * GraphGTest.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter
 * (marvin.ritter@gmail.com)
 */

#include <algorithm>
#include <tuple>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphBuilder.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class GraphGTest : public testing::TestWithParam<std::tuple<bool, bool>> {
public:
    virtual void SetUp();

protected:
    Graph Ghouse;
    std::vector<std::pair<node, node>> houseEdgesOut;
    std::vector<std::vector<edgeweight>> Ahouse;
    count n_house;
    count m_house;

    bool isGraph() const { return !isWeighted() && !isDirected(); }
    bool isWeightedGraph() const { return isWeighted() && !isDirected(); }
    bool isDirectedGraph() const { return !isWeighted() && isDirected(); }
    bool isWeightedDirectedGraph() const { return isWeighted() && isDirected(); }

    bool isWeighted() const;
    bool isDirected() const;
    Graph createGraph(count n = 0) const;
    Graph createGraph(count n, count m) const;
    count countSelfLoopsManually(const Graph &G);
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, GraphGTest, testing::Values(std::make_tuple(false, false),
                        std::make_tuple(true, false),
                        std::make_tuple(false, true),
                        std::make_tuple(true, true)));

bool GraphGTest::isWeighted() const { return std::get<0>(GetParam()); }
bool GraphGTest::isDirected() const { return std::get<1>(GetParam()); }

Graph GraphGTest::createGraph(count n) const {
    bool weighted, directed;
    std::tie(weighted, directed) = GetParam();
    Graph G(n, weighted, directed);
    return G;
}

Graph GraphGTest::createGraph(count n, count m) const {
    auto G = createGraph(n);
    while(G.numberOfEdges() < m) {
        const auto u = Aux::Random::index(n);
        const auto v = Aux::Random::index(n);
        if (u == v) continue;
        if (G.hasEdge(u, v)) continue;

        const auto p = Aux::Random::probability();
        G.addEdge(u, v, p);
    }
    return G;
}

count GraphGTest::countSelfLoopsManually(const Graph &G) {
    count c = 0;
    G.forEdges([&](node u, node v) {
        if (u == v) {
            c += 1;
        }
    });
    return c;
}

void GraphGTest::SetUp() {
    /*
     *    0
     *   . \
     *  /   \
     * /     .
     * 1 <-- 2
     * ^ \  .|
     * |  \/ |
     * | / \ |
     * |/   ..
     * 3 <-- 4
     *
     * move you pen from node to node:
     * 3 -> 1 -> 0 -> 2 -> 1 -> 4 -> 3 -> 2 -> 4
     */
    n_house = 5;
    m_house = 8;

    Ghouse = createGraph(5);
    houseEdgesOut = {{3, 1}, {1, 0}, {0, 2}, {2, 1},
                     {1, 4}, {4, 3}, {3, 2}, {2, 4}};
    Ahouse = {n_house, std::vector<edgeweight>(n_house, 0.0)};
    edgeweight ew = 1.0;
    for (auto &e : houseEdgesOut) {
        node u = e.first;
        node v = e.second;
        Ghouse.addEdge(u, v, ew);

        Ahouse[u][v] = ew;

        if (!Ghouse.isDirected()) {
            Ahouse[v][u] = ew;
        }

        if (Ghouse.isWeighted()) {
            ew += 1.0;
        }
    }
}

/** CONSTRUCTORS **/

TEST_P(GraphGTest, testCopyConstructor) {
    Graph G = Graph(this->Ghouse, false, false);
    Graph GW = Graph(this->Ghouse, true, false);
    Graph D = Graph(this->Ghouse, false, true);
    Graph DW = Graph(this->Ghouse, true, true);

    ASSERT_FALSE(G.isWeighted());
    ASSERT_FALSE(G.isDirected());
    ASSERT_EQ(this->Ghouse.numberOfNodes(), G.numberOfNodes());
    ASSERT_EQ(this->Ghouse.numberOfEdges(), G.numberOfEdges());

    ASSERT_TRUE(GW.isWeighted());
    ASSERT_FALSE(GW.isDirected());
    ASSERT_EQ(this->Ghouse.numberOfNodes(), GW.numberOfNodes());
    ASSERT_EQ(this->Ghouse.numberOfEdges(), GW.numberOfEdges());

    ASSERT_FALSE(D.isWeighted());
    ASSERT_TRUE(D.isDirected());
    ASSERT_EQ(this->Ghouse.numberOfNodes(), D.numberOfNodes());
    ASSERT_EQ(this->Ghouse.numberOfEdges(), D.numberOfEdges());

    ASSERT_TRUE(DW.isWeighted());
    ASSERT_TRUE(DW.isDirected());
    ASSERT_EQ(this->Ghouse.numberOfNodes(), DW.numberOfNodes());
    ASSERT_EQ(this->Ghouse.numberOfEdges(), DW.numberOfEdges());

    this->Ghouse.forNodes([&](node v) {
        count d = this->Ghouse.degree(v);
        count dUndirected = isDirected() ? d + this->Ghouse.degreeIn(v) : d;
        ASSERT_EQ(dUndirected, G.degree(v));
        ASSERT_EQ(dUndirected, GW.degree(v));
        ASSERT_EQ(d, D.degree(v));
        ASSERT_EQ(d, DW.degree(v));
    });

    // if Ghouse was directed we should have an exact copy of it, but if it was
    // undirected we should have edges in both directions
    count m = 0;
    G.forEdges([&](node u, node v) {
        ASSERT_TRUE(G.hasEdge(v, u));
        ASSERT_EQ(defaultEdgeWeight, G.weight(v, u));
        ASSERT_EQ(defaultEdgeWeight, G.weight(u, v));

        auto e = std::make_pair(u, v);
        bool found =
            std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) !=
            this->houseEdgesOut.end();
        if (!found) {
            e = std::make_pair(v, u);
            found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
                              e) != this->houseEdgesOut.end();
        }
        ASSERT_TRUE(found);
        m++;
    });
    ASSERT_EQ(8u, m);

    m = 0;
    GW.forEdges([&](node u, node v) {
        ASSERT_TRUE(GW.hasEdge(v, u));
        ASSERT_EQ(GW.weight(u, v), GW.weight(v, u));

        auto e = std::make_pair(u, v);
        bool found =
            std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) !=
            this->houseEdgesOut.end();
        if (!found) {
            e = std::make_pair(v, u);
            found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
                              e) != this->houseEdgesOut.end();
            ASSERT_EQ(this->Ghouse.weight(v, u), GW.weight(v, u));
        } else {
            ASSERT_EQ(this->Ghouse.weight(u, v), GW.weight(u, v));
        }
        ASSERT_TRUE(found);
        m++;
    });
    ASSERT_EQ(8u, m);

    m = 0;
    D.forEdges([&](node u, node v) {
        ASSERT_EQ(defaultEdgeWeight, D.weight(u, v));

        auto e = std::make_pair(u, v);
        bool found =
            std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) !=
            this->houseEdgesOut.end();
        if (!this->Ghouse.isDirected()) {
            ASSERT_TRUE(D.hasEdge(v, u));

            e = std::make_pair(v, u);
            found = found ||
                    (std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
                               e) != this->houseEdgesOut.end());
        } else {
            ASSERT_FALSE(D.hasEdge(v, u));
        }
        ASSERT_TRUE(found);
        m++;
    });
    count m_expected = isDirected() ? 8 : 16;
    ASSERT_EQ(m_expected, m);

    m = 0;
    DW.forEdges([&](node u, node v) {
        ASSERT_EQ(this->Ghouse.weight(u, v), DW.weight(u, v));

        auto e = std::make_pair(u, v);
        bool found =
            std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) !=
            this->houseEdgesOut.end();
        if (!this->Ghouse.isDirected()) {
            ASSERT_TRUE(DW.hasEdge(v, u));
            e = std::make_pair(v, u);
            found = found ||
                    (std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
                               e) != this->houseEdgesOut.end());
        } else {
            ASSERT_FALSE(DW.hasEdge(v, u));
        }
        ASSERT_TRUE(found);
        m++;
    });
    m_expected = isDirected() ? 8 : 16;
    ASSERT_EQ(m_expected, m);
}

/** NODE MODIFIERS **/

TEST_P(GraphGTest, testAddNode) {
    Graph G = createGraph();

    ASSERT_FALSE(G.hasNode(0));
    ASSERT_FALSE(G.hasNode(1));
    ASSERT_EQ(0u, G.numberOfNodes());

    G.addNode();
    ASSERT_TRUE(G.hasNode(0));
    ASSERT_FALSE(G.hasNode(1));
    ASSERT_EQ(1u, G.numberOfNodes());

    Graph G2 = createGraph(2);
    ASSERT_TRUE(G2.hasNode(0));
    ASSERT_TRUE(G2.hasNode(1));
    ASSERT_FALSE(G2.hasNode(2));
    ASSERT_EQ(2u, G2.numberOfNodes());

    G2.addNode();
    G2.addNode();
    ASSERT_TRUE(G2.hasNode(2));
    ASSERT_TRUE(G2.hasNode(3));
    ASSERT_FALSE(G2.hasNode(4));
    ASSERT_EQ(4u, G2.numberOfNodes());
}

TEST_P(GraphGTest, testAddNodes) {
    auto G = createGraph(5);
    G.addNodes(5);

    ASSERT_EQ(G.numberOfNodes(), 10);
    ASSERT_TRUE(G.hasNode(9));

    G.addNodes(90);

    ASSERT_EQ(G.numberOfNodes(), 100);
    ASSERT_TRUE(G.hasNode(99));
}

TEST_P(GraphGTest, testRemoveNode) {
    auto testGraph = [&](Graph &G) {
        count n = G.numberOfNodes();
        count z = n;
        count m = G.numberOfEdges();
        for (node u = 0; u < z; ++u) {
            count deg = G.degreeOut(u);
            count degIn = G.isDirected() ? G.degreeIn(u) : 0;
            G.removeNode(u);
            --n;
            m -= deg + degIn;
            ASSERT_EQ(G.numberOfNodes(), n);
            ASSERT_EQ(G.numberOfEdges(), m);
            G.forNodes([&](node v) { ASSERT_EQ(G.hasNode(v), v > u); });
        }
    };

    Graph G1 = ErdosRenyiGenerator(200, 0.2, false).generate();
    Graph G2 = ErdosRenyiGenerator(200, 0.2, true).generate();
    testGraph(G1);
    testGraph(G2);
}

TEST_P(GraphGTest, testHasNode) {
    Graph G = createGraph(5);

    ASSERT_TRUE(G.hasNode(0));
    ASSERT_TRUE(G.hasNode(1));
    ASSERT_TRUE(G.hasNode(2));
    ASSERT_TRUE(G.hasNode(3));
    ASSERT_TRUE(G.hasNode(4));
    ASSERT_FALSE(G.hasNode(5));
    ASSERT_FALSE(G.hasNode(6));

    G.removeNode(0);
    G.removeNode(2);
    G.addNode();

    ASSERT_FALSE(G.hasNode(0));
    ASSERT_TRUE(G.hasNode(1));
    ASSERT_FALSE(G.hasNode(2));
    ASSERT_TRUE(G.hasNode(3));
    ASSERT_TRUE(G.hasNode(4));
    ASSERT_TRUE(G.hasNode(5));
    ASSERT_FALSE(G.hasNode(6));
}

TEST_P(GraphGTest, testRestoreNode) {
    Graph G = createGraph(4);

    ASSERT_EQ(4u, G.numberOfNodes());
    ASSERT_TRUE(G.hasNode(0));
    ASSERT_TRUE(G.hasNode(1));
    ASSERT_TRUE(G.hasNode(2));
    ASSERT_TRUE(G.hasNode(3));

    G.removeNode(0);

    ASSERT_EQ(3u, G.numberOfNodes());
    ASSERT_FALSE(G.hasNode(0));
    ASSERT_TRUE(G.hasNode(1));
    ASSERT_TRUE(G.hasNode(2));
    ASSERT_TRUE(G.hasNode(3));

    G.restoreNode(0);

    ASSERT_EQ(4u, G.numberOfNodes());
    ASSERT_TRUE(G.hasNode(0));
    ASSERT_TRUE(G.hasNode(1));
    ASSERT_TRUE(G.hasNode(2));
    ASSERT_TRUE(G.hasNode(3));
}

/** NODE PROPERTIES **/

TEST_P(GraphGTest, testDegree) {
    if (isDirected()) {
        ASSERT_EQ(1u, this->Ghouse.degree(0));
        ASSERT_EQ(2u, this->Ghouse.degree(1));
        ASSERT_EQ(2u, this->Ghouse.degree(2));
        ASSERT_EQ(2u, this->Ghouse.degree(3));
        ASSERT_EQ(1u, this->Ghouse.degree(4));
    } else {
        ASSERT_EQ(2u, this->Ghouse.degree(0));
        ASSERT_EQ(4u, this->Ghouse.degree(1));
        ASSERT_EQ(4u, this->Ghouse.degree(2));
        ASSERT_EQ(3u, this->Ghouse.degree(3));
        ASSERT_EQ(3u, this->Ghouse.degree(4));
    }
}

TEST_P(GraphGTest, testDegreeIn) {
    if (isDirected()) {
        ASSERT_EQ(1u, this->Ghouse.degreeIn(0));
        ASSERT_EQ(2u, this->Ghouse.degreeIn(1));
        ASSERT_EQ(2u, this->Ghouse.degreeIn(2));
        ASSERT_EQ(1u, this->Ghouse.degreeIn(3));
        ASSERT_EQ(2u, this->Ghouse.degreeIn(4));
    } else {
        ASSERT_EQ(2u, this->Ghouse.degreeIn(0));
        ASSERT_EQ(4u, this->Ghouse.degreeIn(1));
        ASSERT_EQ(4u, this->Ghouse.degreeIn(2));
        ASSERT_EQ(3u, this->Ghouse.degreeIn(3));
        ASSERT_EQ(3u, this->Ghouse.degreeIn(4));
    }
}

TEST_P(GraphGTest, testDegreeOut) {
    if (isDirected()) {
        ASSERT_EQ(1u, this->Ghouse.degreeOut(0));
        ASSERT_EQ(2u, this->Ghouse.degreeOut(1));
        ASSERT_EQ(2u, this->Ghouse.degreeOut(2));
        ASSERT_EQ(2u, this->Ghouse.degreeOut(3));
        ASSERT_EQ(1u, this->Ghouse.degreeOut(4));
    } else {
        ASSERT_EQ(2u, this->Ghouse.degreeOut(0));
        ASSERT_EQ(4u, this->Ghouse.degreeOut(1));
        ASSERT_EQ(4u, this->Ghouse.degreeOut(2));
        ASSERT_EQ(3u, this->Ghouse.degreeOut(3));
        ASSERT_EQ(3u, this->Ghouse.degreeOut(4));
    }
}

TEST_P(GraphGTest, testIsIsolated) {
    ASSERT_FALSE(this->Ghouse.isIsolated(0));
    ASSERT_FALSE(this->Ghouse.isIsolated(1));
    ASSERT_FALSE(this->Ghouse.isIsolated(2));
    ASSERT_FALSE(this->Ghouse.isIsolated(3));
    ASSERT_FALSE(this->Ghouse.isIsolated(4));

    this->Ghouse.addNode();
    ASSERT_TRUE(this->Ghouse.isIsolated(5));

    this->Ghouse.removeEdge(1, 0);
    ASSERT_FALSE(this->Ghouse.isIsolated(0));

    this->Ghouse.removeEdge(0, 2);
    ASSERT_TRUE(this->Ghouse.isIsolated(0));

    this->Ghouse.addEdge(1, 0);
    ASSERT_FALSE(this->Ghouse.isIsolated(0));
}

TEST_P(GraphGTest, testWeightedDegree) {
    // add self-loop
    this->Ghouse.addEdge(2, 2, 0.75);

    if (isGraph()) {
        ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(0));
        ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.weightedDegree(1));
        ASSERT_EQ(5 * defaultEdgeWeight, this->Ghouse.weightedDegree(2));
        ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(3));
        ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(4));
    }

    if (isWeightedGraph()) {
        ASSERT_EQ(5.0, this->Ghouse.weightedDegree(0));
        ASSERT_EQ(12.0, this->Ghouse.weightedDegree(1));
        ASSERT_EQ(22.75, this->Ghouse.weightedDegree(2));
        ASSERT_EQ(14.0, this->Ghouse.weightedDegree(3));
        ASSERT_EQ(19.0, this->Ghouse.weightedDegree(4));
    }

    if (isDirectedGraph()) {
        // only count outgoing edges
        ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(0));
        ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(1));
        ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(2));
        ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(3));
        ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(4));
    }

    if (isWeightedDirectedGraph()) {
        // only sum weight of outgoing edges
        ASSERT_EQ(3.0, this->Ghouse.weightedDegree(0));
        ASSERT_EQ(7.0, this->Ghouse.weightedDegree(1));
        ASSERT_EQ(12.75, this->Ghouse.weightedDegree(2));
        ASSERT_EQ(8.0, this->Ghouse.weightedDegree(3));
        ASSERT_EQ(6.0, this->Ghouse.weightedDegree(4));
    }
}

TEST_P(GraphGTest, testWeightedDegree2) {
    // add self-loop
    this->Ghouse.addEdge(2, 2, 0.75);

    if (isGraph()) {
        ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(0, true));
        ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.weightedDegree(1, true));
        ASSERT_EQ(6 * defaultEdgeWeight, this->Ghouse.weightedDegree(2, true));
        ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(3, true));
        ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(4, true));
    }

    if (isWeightedGraph()) {
        ASSERT_EQ(5.0, this->Ghouse.weightedDegree(0, true));
        ASSERT_EQ(12.0, this->Ghouse.weightedDegree(1, true));
        ASSERT_EQ(23.5, this->Ghouse.weightedDegree(2, true));
        ASSERT_EQ(14.0, this->Ghouse.weightedDegree(3, true));
        ASSERT_EQ(19.0, this->Ghouse.weightedDegree(4, true));
    }

    if (isDirectedGraph()) {
        // only count outgoing edges
        ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(0, true));
        ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(1, true));
        ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.weightedDegree(2, true));
        ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(3, true));
        ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(4, true));
    }

    if (isWeightedDirectedGraph()) {
        // only sum weight of outgoing edges
        ASSERT_EQ(3.0, this->Ghouse.weightedDegree(0, true));
        ASSERT_EQ(7.0, this->Ghouse.weightedDegree(1, true));
        ASSERT_EQ(13.5, this->Ghouse.weightedDegree(2, true));
        ASSERT_EQ(8.0, this->Ghouse.weightedDegree(3, true));
        ASSERT_EQ(6.0, this->Ghouse.weightedDegree(4, true));
    }
}

TEST_P(GraphGTest, testWeightedDegree3) {
    constexpr count n = 100;
    constexpr double p = 0.1;

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(n, p, isDirected()).generate();
        if (isWeighted()) {
            G = Graph(G, true, G.isDirected());
            G.forEdges([&](node u, node v) { G.setWeight(u, v, Aux::Random::probability()); });
        }
        G.forNodes([&](node u) {
            edgeweight wDeg = 0, wDegTwice = 0;
            G.forNeighborsOf(u, [&](node v, edgeweight w) {
                wDeg += w;
                wDegTwice += (u == v) ? 2. * w : w;
            });

            EXPECT_DOUBLE_EQ(G.weightedDegree(u), wDeg);
            EXPECT_DOUBLE_EQ(G.weightedDegree(u, true), wDegTwice);

            edgeweight wInDeg = 0, wInDegTwice = 0;
            G.forInNeighborsOf(u, [&](node v, edgeweight w) {
                wInDeg += w;
                wInDegTwice += (u == v) ? 2. * w : w;
            });

            EXPECT_DOUBLE_EQ(G.weightedDegreeIn(u), wInDeg);
            EXPECT_DOUBLE_EQ(G.weightedDegreeIn(u, true), wInDegTwice);
        });
    }
}

/** EDGE MODIFIERS **/

TEST_P(GraphGTest, testAddEdge) {
    Graph G = createGraph(3);

    // Graph without edges
    ASSERT_EQ(0u, G.numberOfEdges());
    ASSERT_FALSE(G.hasEdge(0, 2));
    ASSERT_FALSE(G.hasEdge(0, 1));
    ASSERT_FALSE(G.hasEdge(1, 2));
    ASSERT_FALSE(G.hasEdge(2, 2));
    ASSERT_EQ(nullWeight, G.weight(0, 2));
    ASSERT_EQ(nullWeight, G.weight(0, 1));
    ASSERT_EQ(nullWeight, G.weight(1, 2));
    ASSERT_EQ(nullWeight, G.weight(2, 2));

    // Graph with 2 normal edges
    G.addEdge(0, 1, 4.51);
    G.addEdge(1, 2, 2.39);
    ASSERT_EQ(2u, G.numberOfEdges());
    ASSERT_FALSE(G.hasEdge(0, 2)); // was never added
    ASSERT_TRUE(G.hasEdge(0, 1));
    ASSERT_TRUE(G.hasEdge(1, 2));
    ASSERT_FALSE(G.hasEdge(2, 2)); // will be added later

    // check weights
    if (G.isWeighted()) {
        ASSERT_EQ(4.51, G.weight(0, 1));
        ASSERT_EQ(2.39, G.weight(1, 2));
    } else {
        ASSERT_EQ(defaultEdgeWeight, G.weight(0, 1));
        ASSERT_EQ(defaultEdgeWeight, G.weight(1, 2));
    }

    if (G.isDirected()) {
        ASSERT_FALSE(G.hasEdge(1, 0));
        ASSERT_FALSE(G.hasEdge(2, 1));

        // add edge in the other direction
        // note: bidirectional edges are not supported, so both edges have different
        // weights
        G.addEdge(2, 1, 6.23);
        ASSERT_TRUE(G.hasEdge(2, 1));
        if (G.isWeighted()) {
            ASSERT_EQ(2.39, G.weight(1, 2));
            ASSERT_EQ(6.23, G.weight(2, 1));
        } else {
            ASSERT_EQ(defaultEdgeWeight, G.weight(2, 1));
        }
    } else {
        ASSERT_TRUE(G.hasEdge(1, 0));
        ASSERT_TRUE(G.hasEdge(2, 1));
        if (G.isWeighted()) {
            ASSERT_EQ(4.51, G.weight(1, 0));
            ASSERT_EQ(2.39, G.weight(2, 1));
        } else {
            ASSERT_EQ(defaultEdgeWeight, G.weight(1, 0));
            ASSERT_EQ(defaultEdgeWeight, G.weight(2, 1));
        }
    }

    // add self loop
    G.addEdge(2, 2, 0.72);
    ASSERT_TRUE(G.hasEdge(2, 2));
    if (G.isWeighted()) {
        ASSERT_EQ(0.72, G.weight(2, 2));
    } else {
        ASSERT_EQ(defaultEdgeWeight, G.weight(2, 2));
    }
}

TEST_P(GraphGTest, testRemoveEdge) {
    double epsilon = 1e-6;
    Graph G = createGraph(3);

    edgeweight ewBefore = G.totalEdgeWeight();

    G.addEdge(0, 1, 3.14);

    if (G.isWeighted()) {
        ASSERT_NEAR(ewBefore + 3.14, G.totalEdgeWeight(), epsilon);
    } else {
        ASSERT_NEAR(ewBefore + defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
    }

    G.addEdge(0, 0);

    ASSERT_EQ(2u, G.numberOfEdges());
    ASSERT_TRUE(G.hasEdge(0, 0));
    ASSERT_TRUE(G.hasEdge(0, 1));
    ASSERT_FALSE(G.hasEdge(2, 1));

    // test remove regular edge
    ewBefore = G.totalEdgeWeight();
    G.removeEdge(0, 1);
    if (G.isWeighted()) {
        ASSERT_NEAR(ewBefore - 3.14, G.totalEdgeWeight(), epsilon);
    } else {
        ASSERT_NEAR(ewBefore - defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
    }

    ASSERT_EQ(1u, G.numberOfEdges());
    ASSERT_TRUE(G.hasEdge(0, 0));
    ASSERT_FALSE(G.hasEdge(0, 1));
    ASSERT_FALSE(G.hasEdge(2, 1));

    // test remove self-loop
    G.addEdge(2, 1);

    ewBefore = G.totalEdgeWeight();
    G.removeEdge(0, 0);
    if (G.isWeighted()) {
        ASSERT_NEAR(ewBefore - defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
    } else {
        ASSERT_NEAR(ewBefore - defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
    }

    ASSERT_EQ(1u, G.numberOfEdges());
    ASSERT_FALSE(G.hasEdge(0, 0));
    ASSERT_FALSE(G.hasEdge(0, 1));
    ASSERT_TRUE(G.hasEdge(2, 1));

    // test from removeselfloops adapted for removeEdge
    G = createGraph(2);

    ewBefore = G.totalEdgeWeight();

    G.addEdge(0, 1);
    G.addEdge(0, 0, 3.14);
    G.addEdge(1, 1);

    if (G.isWeighted()) {
        EXPECT_NEAR(ewBefore + 3.14 + 2 * defaultEdgeWeight, G.totalEdgeWeight(),
                    epsilon);
    } else {
        EXPECT_NEAR(ewBefore + 3 * defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
    }

    EXPECT_EQ(3u, G.numberOfEdges());
    EXPECT_TRUE(G.hasEdge(0, 0));
    EXPECT_TRUE(G.hasEdge(0, 1));
    EXPECT_TRUE(G.hasEdge(1, 1));
    EXPECT_EQ(G.numberOfSelfLoops(), 2u);

    // remove self-loops
    ewBefore = G.totalEdgeWeight();

    G.removeEdge(0, 0);
    G.removeEdge(1, 1);

    if (G.isWeighted()) {
        EXPECT_NEAR(ewBefore - defaultEdgeWeight - 3.14, G.totalEdgeWeight(),
                    epsilon);
    } else {
        EXPECT_NEAR(ewBefore - 2 * defaultEdgeWeight, G.totalEdgeWeight(), epsilon)
            << "Weighted, directed: " << G.isWeighted() << ", " << G.isDirected();
    }

    EXPECT_EQ(1u, G.numberOfEdges());
    EXPECT_FALSE(G.hasEdge(0, 0));
    EXPECT_FALSE(G.hasEdge(1, 1));
    EXPECT_TRUE(G.hasEdge(0, 1));
    EXPECT_EQ(0u, G.numberOfSelfLoops())
        << "Weighted, directed: " << G.isWeighted() << ", " << G.isDirected();
}

TEST_P(GraphGTest, testRemoveAllEdges) {
    constexpr count n = 100;
    constexpr double p = 0.2;

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto g = ErdosRenyiGenerator(n, p, isDirected()).generate();
        if (isWeighted()) {
            g = Graph(g, true, isDirected());
        }

        g.removeAllEdges();

        EXPECT_EQ(g.numberOfEdges(), 0);

        count edgeCount = 0;
        g.forEdges([&edgeCount](node, node) { ++edgeCount; });
        EXPECT_EQ(edgeCount, 0);

        g.forNodes([&](node u) {
            EXPECT_EQ(g.degree(u), 0);
            EXPECT_EQ(g.degree(u), 0);
        });
    }
}

TEST_P(GraphGTest, testRemoveSelfLoops) {
    constexpr count n = 100;
    constexpr count nSelfLoops = 100;
    constexpr double p = 0.2;

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto g = ErdosRenyiGenerator(n, p, isDirected()).generate();
        if (isWeighted()) {
            g = Graph(g, true, isDirected());
        }

        for (count i = 0; i < nSelfLoops; ++i) {
            const auto u = GraphTools::randomNode(g);
            g.addEdge(u, u);
        }

        const auto numberOfSelfLoops = g.numberOfSelfLoops();
        const auto numberOfEdges = g.numberOfEdges();
        g.removeSelfLoops();

        EXPECT_EQ(numberOfEdges - numberOfSelfLoops, g.numberOfEdges());
        EXPECT_EQ(g.numberOfSelfLoops(), 0);
        g.forNodes([&g](const node u) { EXPECT_FALSE(g.hasEdge(u, u)); });
    }
}

TEST_P(GraphGTest, testRemoveMultiEdges) {
    constexpr count n = 200;
    constexpr double p = 0.1;
    constexpr count nMultiEdges = 10;
    constexpr count nMultiSelfLoops = 10;

    auto getGraphEdges = [](const Graph &G) {
        std::vector<std::pair<node, node>> edges;
        edges.reserve(G.numberOfEdges());

        G.forEdges([&](const node u, const node v) {
            edges.push_back({u, v});
        });

        return edges;
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto g = ErdosRenyiGenerator(n, p, isDirected()).generate();
        if (isWeighted()) {
            g = Graph(g, true, isDirected());
        }

        const auto edgeSet = getGraphEdges(g);
        const auto m = g.numberOfEdges();

        // Adding multiedges at random
        for (count i = 0; i < nMultiEdges; ++i) {
            const auto e = GraphTools::randomEdge(g);
            g.addEdge(e.first, e.second);
        }

        std::unordered_set<node> uniqueSelfLoops;
        // Adding multiple self-loops at random
        for (count i = 0; i < nMultiSelfLoops; ++i) {
            const auto u = GraphTools::randomNode(g);
            g.addEdge(u, u);
            g.addEdge(u, u);
            uniqueSelfLoops.insert(u);
        }

        EXPECT_EQ(g.numberOfEdges(), m + nMultiEdges + 2 * nMultiSelfLoops);

        g.removeMultiEdges();

        EXPECT_EQ(g.numberOfEdges(), m + uniqueSelfLoops.size());
        g.removeSelfLoops();

        EXPECT_EQ(g.numberOfEdges(), m);
        auto edgeSet_ = getGraphEdges(g);

        for (count i = 0; i < g.numberOfEdges(); ++i)
            EXPECT_EQ(edgeSet[i], edgeSet_[i]);
    }
}

TEST_P(GraphGTest, testHasEdge) {
    auto containsEdge = [&](std::pair<node, node> e) {
        auto it =
            std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
        return it != this->houseEdgesOut.end();
    };

    for (node u = 0; u < this->Ghouse.upperNodeIdBound(); u++) {
        for (node v = 0; v < this->Ghouse.upperNodeIdBound(); v++) {
            auto edge = std::make_pair(u, v);
            auto edgeReverse = std::make_pair(v, u);
            bool hasEdge = containsEdge(edge);
            bool hasEdgeReverse = containsEdge(edgeReverse);
            if (this->Ghouse.isDirected()) {
                ASSERT_EQ(hasEdge, this->Ghouse.hasEdge(u, v));
            } else {
                ASSERT_EQ(hasEdge || hasEdgeReverse, this->Ghouse.hasEdge(u, v));
            }
        }
    }
}

/** GLOBAL PROPERTIES **/

TEST_P(GraphGTest, testSelfLoopCountSimple) {
    Graph G(Ghouse);
    G.addEdge(0, 0);
    EXPECT_EQ(1, G.numberOfSelfLoops());
}

TEST_P(GraphGTest, testIsWeighted) {
    ASSERT_EQ(isWeighted(), this->Ghouse.isWeighted());
}

TEST_P(GraphGTest, testIsDirected) {
    ASSERT_EQ(isDirected(), this->Ghouse.isDirected());
}

TEST_P(GraphGTest, testIsEmpty) {
    Graph G1 = createGraph(0);
    Graph G2 = createGraph(2);

    ASSERT_TRUE(G1.isEmpty());
    ASSERT_FALSE(G2.isEmpty());

    node v = G1.addNode();
    G2.removeNode(GraphTools::randomNode(G2));
    ASSERT_FALSE(G1.isEmpty());
    ASSERT_FALSE(G2.isEmpty());

    G1.removeNode(v);
    G2.removeNode(GraphTools::randomNode(G2));
    ASSERT_TRUE(G1.isEmpty());
    ASSERT_TRUE(G2.isEmpty());
}

TEST_P(GraphGTest, testNumberOfNodes) {
    ASSERT_EQ(this->n_house, this->Ghouse.numberOfNodes());

    Graph G1 = createGraph(0);
    ASSERT_EQ(0u, G1.numberOfNodes());
    G1.addNode();
    ASSERT_EQ(1u, G1.numberOfNodes());
    G1.addNode();
    ASSERT_EQ(2u, G1.numberOfNodes());
    G1.removeNode(0);
    ASSERT_EQ(1u, G1.numberOfNodes());
    G1.removeNode(1);
    ASSERT_EQ(0u, G1.numberOfNodes());
}

TEST_P(GraphGTest, testNumberOfEdges) {
    ASSERT_EQ(this->m_house, this->Ghouse.numberOfEdges());

    Graph G1 = createGraph(5);
    ASSERT_EQ(0u, G1.numberOfEdges());
    G1.addEdge(0, 1);
    ASSERT_EQ(1u, G1.numberOfEdges());
    G1.addEdge(1, 2);
    ASSERT_EQ(2u, G1.numberOfEdges());
    G1.removeEdge(0, 1);
    ASSERT_EQ(1u, G1.numberOfEdges());
    G1.removeEdge(1, 2);
    ASSERT_EQ(0u, G1.numberOfEdges());
}

TEST_P(GraphGTest, testNumberOfSelfLoops) {
    Graph G = createGraph(3);
    G.addEdge(0, 1);
    ASSERT_EQ(0u, G.numberOfSelfLoops());
    G.addEdge(0, 0);
    ASSERT_EQ(1u, G.numberOfSelfLoops());
    G.addEdge(1, 1);
    G.addEdge(1, 2);
    ASSERT_EQ(2u, G.numberOfSelfLoops());
    G.removeEdge(0, 0);
    ASSERT_EQ(1u, G.numberOfSelfLoops());
}

TEST_P(GraphGTest, testSelfLoopConversion) {
    Aux::Random::setSeed(1, false);
    const count runs = 100;
    const count n_max = 200;
    for (index i = 0; i < runs; i++) {
        bool directed = Aux::Random::probability() < 0.5;
        count n = Aux::Random::integer(n_max);
        Graph G(n, false, directed);

        G.forNodes([&](node v) {
            double p = Aux::Random::probability();

            if (p < 0.1) { // new node
                n++;
                G.addNode();
            } else {                                   // new edge
                node u = Aux::Random::integer(v, n - 1); // self-loops possible
                G.addEdge(v, u);
            }
        });
        count measuredSelfLoops = countSelfLoopsManually(G);
        EXPECT_EQ(G.numberOfSelfLoops(), measuredSelfLoops);
        Graph G_converted(G, false, !directed);
        EXPECT_EQ(G_converted.numberOfSelfLoops(), measuredSelfLoops);
    }
}

TEST_P(GraphGTest, testUpperNodeIdBound) {
    ASSERT_EQ(5u, this->Ghouse.upperNodeIdBound());

    Graph G1 = createGraph(0);
    ASSERT_EQ(0u, G1.upperNodeIdBound());
    G1.addNode();
    ASSERT_EQ(1u, G1.upperNodeIdBound());
    G1.addNode();
    ASSERT_EQ(2u, G1.upperNodeIdBound());
    G1.removeNode(1);
    ASSERT_EQ(2u, G1.upperNodeIdBound());
    G1.addNode();
    ASSERT_EQ(3u, G1.upperNodeIdBound());
}

TEST_P(GraphGTest, testCheckConsistency_MultiEdgeDetection) {
    Graph G = createGraph(3);
    ASSERT_TRUE(G.checkConsistency());
    G.addEdge(0, 1);
    ASSERT_TRUE(G.checkConsistency());
    G.addEdge(0, 2);
    G.addEdge(0, 1);
    ASSERT_FALSE(G.checkConsistency());
    G.removeEdge(0, 1);
    ASSERT_TRUE(G.checkConsistency());
    G.removeEdge(0, 1);
    ASSERT_TRUE(G.checkConsistency());
}

/** EDGE ATTRIBUTES **/

TEST_P(GraphGTest, testWeight) {
    this->Ghouse.forNodes([&](node u) {
        this->Ghouse.forNodes([&](node v) {
            ASSERT_EQ(this->Ahouse[u][v], this->Ghouse.weight(u, v));
        });
    });
}

TEST_P(GraphGTest, testSetWeight) {
    Graph G = createGraph(10);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    if (isWeighted()) {
        // edges should get weight defaultWeight on creation and setWeight should
        // overwrite this
        G.setWeight(1, 2, 2.718);
        EXPECT_EQ(defaultEdgeWeight, G.weight(0, 1));
        EXPECT_EQ(2.718, G.weight(1, 2));
        if (isDirected()) {
            EXPECT_EQ(nullWeight, G.weight(1, 0));
            EXPECT_EQ(nullWeight, G.weight(2, 1));
        } else {
            // undirected graph is symmetric
            EXPECT_EQ(defaultEdgeWeight, G.weight(1, 0));
            EXPECT_EQ(2.718, G.weight(2, 1));
        }

        // setting an edge weight should create the edge if it doesn't exists
        ASSERT_FALSE(G.hasEdge(5, 6));
        G.setWeight(5, 6, 56.0);
        ASSERT_EQ(56.0, G.weight(5, 6));
        ASSERT_EQ(isDirected() ? nullWeight : 56.0, G.weight(6, 5));
        ASSERT_TRUE(G.hasEdge(5, 6));

        // directed graphs are not symmetric, undirected are
        G.setWeight(2, 1, 5.243);
        if (isDirected()) {
            EXPECT_EQ(2.718, G.weight(1, 2));
            EXPECT_EQ(5.243, G.weight(2, 1));
        } else {
            EXPECT_EQ(5.243, G.weight(1, 2));
            EXPECT_EQ(5.243, G.weight(2, 1));
        }

        // self-loop
        G.addEdge(4, 4, 2.5);
        ASSERT_EQ(2.5, G.weight(4, 4));
        G.setWeight(4, 4, 3.14);
        ASSERT_EQ(3.14, G.weight(4, 4));
    } else {
        EXPECT_ANY_THROW(G.setWeight(0, 1, 1.5));
    }
}

TEST_P(GraphGTest, increaseWeight) {
    Graph G = createGraph(5);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(3, 4, 3.14);

    if (G.isWeighted()) {
        G.increaseWeight(1, 2, 0.5);
        G.increaseWeight(3, 4, -0.5);

        ASSERT_EQ(defaultEdgeWeight, G.weight(0, 1));
        ASSERT_EQ(defaultEdgeWeight + 0.5, G.weight(1, 2));
        ASSERT_EQ(3.14 - 0.5, G.weight(3, 4));

        if (G.isDirected()) {
            // reverse edges do net exist => weight should be nullWeight
            ASSERT_EQ(nullWeight, G.weight(1, 0));
            ASSERT_EQ(nullWeight, G.weight(2, 1));
            ASSERT_EQ(nullWeight, G.weight(4, 3));
        } else {
            ASSERT_EQ(defaultEdgeWeight, G.weight(1, 0));
            ASSERT_EQ(defaultEdgeWeight + 0.5, G.weight(2, 1));
            ASSERT_EQ(3.14 - 0.5, G.weight(3, 4));
        }
    } else {
        EXPECT_ANY_THROW(G.increaseWeight(1, 2, 0.3));
        EXPECT_ANY_THROW(G.increaseWeight(2, 3, 0.3)); // edge does not exists
    }
}

/** SUMS **/

TEST_P(GraphGTest, testTotalEdgeWeight) {
    Graph G1 = createGraph(5);
    Graph G2 = createGraph(5);
    G2.addEdge(0, 1, 3.14);

    if (this->Ghouse.isWeighted()) {
        ASSERT_EQ(0.0, G1.totalEdgeWeight());
        ASSERT_EQ(3.14, G2.totalEdgeWeight());
        ASSERT_EQ(36.0, this->Ghouse.totalEdgeWeight());
    } else {
        ASSERT_EQ(0 * defaultEdgeWeight, G1.totalEdgeWeight());
        ASSERT_EQ(1 * defaultEdgeWeight, G2.totalEdgeWeight());
        ASSERT_EQ(8 * defaultEdgeWeight, this->Ghouse.totalEdgeWeight());
    }
}

/** Collections **/

TEST_P(GraphGTest, testNodeIterator) {
    Aux::Random::setSeed(42, false);

    auto testForward = [](const Graph &G) {
        auto preIter = G.nodeRange().begin();
        auto postIter = G.nodeRange().begin();

        G.forNodes([&](const node u) {
            ASSERT_EQ(*preIter, u);
            ASSERT_EQ(*postIter, u);
            ++preIter;
            postIter++;
        });

        ASSERT_EQ(preIter, G.nodeRange().end());
        ASSERT_EQ(postIter, G.nodeRange().end());

        Graph G1(G);

        for (const auto u : Graph::NodeRange(G)) {
            ASSERT_TRUE(G1.hasNode(u));
            G1.removeNode(u);
        }

        ASSERT_EQ(G1.numberOfNodes(), 0);
    };

    auto testBackward = [](const Graph &G) {
        const std::vector<node> nodes(Graph::NodeRange(G).begin(), Graph::NodeRange(G).end());
        std::vector<node> v;
        G.forNodes([&](node u) {v.push_back(u);});

        ASSERT_EQ(std::unordered_set<node>(nodes.begin(), nodes.end()).size(), nodes.size());
        ASSERT_EQ(nodes.size(), G.numberOfNodes());

        auto preIter = G.nodeRange().begin();
        auto postIter = G.nodeRange().begin();
        for (count i = 0; i < G.numberOfNodes(); ++i) {
            ++preIter;
            postIter++;
        }

        ASSERT_EQ(preIter, G.nodeRange().end());
        ASSERT_EQ(postIter, G.nodeRange().end());
        auto vecIter = nodes.rbegin();
        while (vecIter != nodes.rend()) {
            ASSERT_EQ(*vecIter, *(--preIter));
            if (postIter != G.nodeRange().end()) {
                ASSERT_NE(*vecIter, *(postIter--));
            } else {
                postIter--;
            }
            ASSERT_EQ(*vecIter, *postIter);
            ++vecIter;
        }

        ASSERT_EQ(preIter, G.nodeRange().begin());
        ASSERT_EQ(postIter, G.nodeRange().begin());
    };

    Graph G(this->Ghouse);
    testForward(G);
    testBackward(G);

    G.removeNode(GraphTools::randomNode(G));
    G.removeNode(GraphTools::randomNode(G));

    testForward(G);
    testBackward(G);
}

TEST_P(GraphGTest, testEdgeIterator) {
    Graph G(this->Ghouse);

    auto testForward = [&](const Graph &G) {
        Graph G1(G);
        auto preIter = G.edgeRange().begin();
        auto postIter = G.edgeRange().begin();

        G.forEdges([&](node, node) {
            ASSERT_EQ(preIter, postIter);
            const auto edge = *preIter;
            ASSERT_TRUE(G.hasEdge(edge.u, edge.v));
            G1.removeEdge(edge.u, edge.v);
            ++preIter;
            postIter++;
        });

        ASSERT_EQ(G1.numberOfEdges(), 0);
        ASSERT_EQ(preIter, G.edgeRange().end());
        ASSERT_EQ(postIter, G.edgeRange().end());

        G1 = G;
        for (const auto edge : Graph::EdgeRange(G)) {
            ASSERT_TRUE(G1.hasEdge(edge.u, edge.v));
            G1.removeEdge(edge.u, edge.v);
        }

        ASSERT_EQ(G1.numberOfEdges(), 0);
    };

    auto testForwardWeighted = [&](const Graph &G) {
        Graph G1(G);
        auto preIter = G.edgeWeightRange().begin();
        auto postIter = preIter;

        G.forEdges([&](node, node) {
            ASSERT_EQ(preIter, postIter);

            const auto edge = *preIter;
            ASSERT_TRUE(G.hasEdge(edge.u, edge.v));
            ASSERT_DOUBLE_EQ(G.weight(edge.u, edge.v), edge.weight);
            G1.removeEdge(edge.u, edge.v);
            ++preIter;
            postIter++;
        });

        ASSERT_EQ(G1.numberOfEdges(), 0);
        ASSERT_EQ(preIter, G.edgeWeightRange().end());
        ASSERT_EQ(postIter, G.edgeWeightRange().end());

        G1 = G;
        for (const auto edge : Graph::EdgeWeightRange(G)) {
            ASSERT_TRUE(G1.hasEdge(edge.u, edge.v));
            ASSERT_DOUBLE_EQ(G1.weight(edge.u, edge.v), edge.weight);
            G1.removeEdge(edge.u, edge.v);
        }

        ASSERT_EQ(G1.numberOfEdges(), 0);
    };

    auto testBackward = [&](const Graph &G) {
        Graph G1(G);
        auto preIter = G.edgeRange().begin();
        auto postIter = preIter;
        G.forEdges([&](node, node) {
            ++preIter;
            postIter++;
        });

        ASSERT_EQ(preIter, G.edgeRange().end());
        ASSERT_EQ(postIter, G.edgeRange().end());

        G.forEdges([&](node, node) {
            --preIter;
            postIter--;
            ASSERT_EQ(preIter, postIter);
            const auto edge = *preIter;
            ASSERT_TRUE(G.hasEdge(edge.u, edge.v));
            G1.removeEdge(edge.u, edge.v);
        });

        ASSERT_EQ(G1.numberOfEdges(), 0);
    };

    auto testBackwardWeighted = [&](const Graph &G) {
        Graph G1(G);
        auto preIter = G.edgeWeightRange().begin();
        auto postIter = preIter;
        G.forEdges([&](node, node) {
            ++preIter;
            postIter++;
        });

        G.forEdges([&](node, node) {
            --preIter;
            postIter--;
            ASSERT_EQ(preIter, postIter);

            const auto edge = *preIter;
            ASSERT_TRUE(G.hasEdge(edge.u, edge.v));
            ASSERT_DOUBLE_EQ(G.weight(edge.u, edge.v), edge.weight);
            G1.removeEdge(edge.u, edge.v);
        });

        ASSERT_EQ(G1.numberOfEdges(), 0);
    };

    auto doTests = [&](const Graph &G) {
        testForward(G);
        testBackward(G);
        testForwardWeighted(G);
        testBackwardWeighted(G);
    };

    doTests(G);

    for (int seed : {1, 2, 3, 4, 5}) {
        Aux::Random::setSeed(seed, false);
        Graph G1(G);
        for (int i = 0; i < 3; ++i) {
            auto e = GraphTools::randomEdge(G1);
            G1.removeEdge(e.first, e.second);
        }

        doTests(G1);
    }
}

TEST_P(GraphGTest, testNeighborsIterators) {
    auto iter = this->Ghouse.neighborRange(1).begin();
    this->Ghouse.forNeighborsOf(1, [&](node v) {
        ASSERT_TRUE(*iter == v);
        ++iter;
    });
    ASSERT_TRUE(iter == this->Ghouse.neighborRange(1).end());

    if (this->Ghouse.isWeighted()) {
        auto iterW = this->Ghouse.weightNeighborRange(1).begin();
        this->Ghouse.forNeighborsOf(1, [&](node v, edgeweight w) {
            ASSERT_TRUE((*iterW).first == v);
            ASSERT_TRUE((*iterW).second == w);
            ++iterW;
        });
        ASSERT_TRUE(iterW == this->Ghouse.weightNeighborRange(1).end());
    }

    if (this->Ghouse.isDirected()) {
        auto inIter = this->Ghouse.inNeighborRange(1).begin();
        this->Ghouse.forInNeighborsOf(1, [&](node v) {
            ASSERT_TRUE(*inIter == v);
            ++inIter;
        });
        ASSERT_TRUE(inIter == this->Ghouse.inNeighborRange(1).end());

        if (this->Ghouse.isWeighted()) {
            auto iterW = this->Ghouse.weightInNeighborRange(1).begin();
            this->Ghouse.forInNeighborsOf(1, [&](node v, edgeweight w) {
                ASSERT_TRUE((*iterW).first == v);
                ASSERT_TRUE((*iterW).second == w);
                ++iterW;
            });
            ASSERT_TRUE(iterW == this->Ghouse.weightInNeighborRange(1).end());
        }
    }
}

/** NODE ITERATORS **/

TEST_P(GraphGTest, testForNodes) {
    Graph G = createGraph(3);
    std::vector<bool> visited(4, false);
    G.forNodes([&](node v) {
        ASSERT_FALSE(visited[v]);
        if (v == 2) {
            G.addNode();
        }
        visited[v] = true;
    });
    for (bool b : visited) {
        ASSERT_TRUE(b);
    }
}

TEST_P(GraphGTest, testParallelForNodes) {
    std::vector<node> visited(Ghouse.upperNodeIdBound());
    this->Ghouse.parallelForNodes([&](node u) {
        visited[u] = u;
    });

    Aux::Parallel::sort(visited.begin(), visited.end());

    ASSERT_EQ(5u, visited.size());
    for (index i = 0; i < this->Ghouse.upperNodeIdBound(); i++) {
        ASSERT_EQ(i, visited[i]);
    }
}

TEST_P(GraphGTest, forNodesWhile) {
    count n = 100;
    Graph G = createGraph(n);
    count stopAfter = 10;
    count nodesSeen = 0;

    G.forNodesWhile([&]() { return nodesSeen < stopAfter; },
                    [&](node) { nodesSeen++; });

    ASSERT_EQ(stopAfter, nodesSeen);
}

TEST_P(GraphGTest, testForNodesInRandomOrder) {
    count n = 1000;
    count samples = 100;
    double maxAbsoluteError = 0.005;
    Graph G = createGraph(n);

    node lastNode = n / 2;
    count greaterLastNode = 0;
    count smallerLastNode = 0;
    std::vector<count> visitCount(n, 0);

    for (count i = 0; i < samples; i++) {
        G.forNodesInRandomOrder([&](node v) {
            if (v > lastNode) {
                greaterLastNode++;
            } else {
                smallerLastNode++;
            }
            visitCount[v]++;
            lastNode = v;
        });
    }

    for (node v = 0; v < n; v++) {
        ASSERT_EQ(samples, visitCount[v]);
    }

    ASSERT_NEAR(0.5, (double)greaterLastNode / n / samples, maxAbsoluteError);
    ASSERT_NEAR(0.5, (double)smallerLastNode / n / samples, maxAbsoluteError);
}

TEST_P(GraphGTest, testForNodePairs) {
    count n = 10;
    count m = n * (n - 1) / 2;
    Graph G = createGraph(n);

    // add all edges
    G.forNodePairs([&](node u, node v) {
        ASSERT_FALSE(G.hasEdge(u, v));
        G.addEdge(u, v);
        ASSERT_TRUE(G.hasEdge(u, v));
    });

    EXPECT_EQ(m, G.numberOfEdges());

    // remove all edges
    G.forNodePairs([&](node u, node v) {
        ASSERT_TRUE(G.hasEdge(u, v));
        G.removeEdge(u, v);
        ASSERT_FALSE(G.hasEdge(u, v));
    });

    EXPECT_EQ(0u, G.numberOfEdges());
}

/** EDGE ITERATORS **/

TEST_P(GraphGTest, testForEdges) {
    Graph G = createGraph(4);
    G.addEdge(0, 1); // 0 * 1 = 0
    G.addEdge(1, 2); // 1 * 2 = 2
    G.addEdge(3, 2); // 3 * 2 = 1 (mod 5)
    G.addEdge(2, 2); // 2 * 2 = 4
    G.addEdge(3, 1); // 3 * 1 = 3

    std::vector<bool> edgesSeen(5, false);

    G.forEdges([&](node u, node v) {
        ASSERT_TRUE(G.hasEdge(u, v));
        index id = (u * v) % 5;
        edgesSeen[id] = true;
    });

    for (auto b : edgesSeen) {
        ASSERT_TRUE(b);
    }
}

TEST_P(GraphGTest, testForWeightedEdges) {
    double epsilon = 1e-6;

    Graph G = createGraph(4);
    G.addEdge(0, 1, 0.1); // 0 * 1 = 0
    G.addEdge(3, 2, 0.2); // 3 * 2 = 1 (mod 5)
    G.addEdge(1, 2, 0.3); // 1 * 2 = 2
    G.addEdge(3, 1, 0.4); // 3 * 1 = 3
    G.addEdge(2, 2, 0.5); // 2 * 2 = 4

    std::vector<bool> edgesSeen(5, false);

    edgeweight weightSum = 0;
    G.forEdges([&](node u, node v, edgeweight ew) {
        ASSERT_TRUE(G.hasEdge(u, v));
        ASSERT_EQ(G.weight(u, v), ew);

        index id = (u * v) % 5;
        edgesSeen[id] = true;
        if (G.isWeighted()) {
            ASSERT_NEAR((id + 1) * 0.1, ew, epsilon);
        } else {
            ASSERT_EQ(defaultEdgeWeight, ew);
        }
        weightSum += ew;
    });

    for (auto b : edgesSeen) {
        ASSERT_TRUE(b);
    }
    if (G.isWeighted()) {
        ASSERT_NEAR(1.5, weightSum, epsilon);
    } else {
        ASSERT_NEAR(5 * defaultEdgeWeight, weightSum, epsilon);
    }
}

TEST_P(GraphGTest, testParallelForWeightedEdges) {
    count n = 4;
    Graph G = createGraph(n);
    G.forNodePairs([&](node u, node v) { G.addEdge(u, v, 1.0); });

    edgeweight weightSum = 0.0;
    G.parallelForEdges([&](node, node, edgeweight ew) {
#pragma omp atomic
        weightSum += ew;
    });

    ASSERT_EQ(6.0, weightSum) << "sum of edge weights should be 6 in every case";
}

TEST_P(GraphGTest, testParallelForEdges) {
    count n = 4;
    Graph G = createGraph(n);
    G.forNodePairs([&](node u, node v) { G.addEdge(u, v); });

    edgeweight weightSum = 0.0;
    G.parallelForEdges([&](node, node) {
#pragma omp atomic
        weightSum += 1;
    });

    ASSERT_EQ(6.0, weightSum) << "sum of edge weights should be 6 in every case";
}

/** NEIGHBORHOOD ITERATORS **/

TEST_P(GraphGTest, testForNeighborsOf) {
    std::vector<node> visited;
    this->Ghouse.forNeighborsOf(3, [&](node u) { visited.push_back(u); });

    Aux::Parallel::sort(visited.begin(), visited.end());

    if (isDirected()) {
        ASSERT_EQ(2u, visited.size());
        ASSERT_EQ(1u, visited[0]);
        ASSERT_EQ(2u, visited[1]);
    } else {
        ASSERT_EQ(3u, visited.size());
        ASSERT_EQ(1u, visited[0]);
        ASSERT_EQ(2u, visited[1]);
        ASSERT_EQ(4u, visited[2]);
    }
}

TEST_P(GraphGTest, testForWeightedNeighborsOf) {
    std::vector<std::pair<node, edgeweight>> visited;
    this->Ghouse.forNeighborsOf(3, [&](node u, edgeweight ew) {
        visited.push_back(std::make_pair(u, ew));
    });

    // should sort after the first element
    Aux::Parallel::sort(visited.begin(), visited.end());

    if (isGraph()) {
        ASSERT_EQ(3u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(4u, visited[2].first);
        ASSERT_EQ(defaultEdgeWeight, visited[0].second);
        ASSERT_EQ(defaultEdgeWeight, visited[1].second);
        ASSERT_EQ(defaultEdgeWeight, visited[2].second);
    }

    if (isWeightedGraph()) {
        ASSERT_EQ(3u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(4u, visited[2].first);
        ASSERT_EQ(1.0, visited[0].second);
        ASSERT_EQ(7.0, visited[1].second);
        ASSERT_EQ(6.0, visited[2].second);
    }

    if (isDirectedGraph()) {
        ASSERT_EQ(2u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(defaultEdgeWeight, visited[0].second);
        ASSERT_EQ(defaultEdgeWeight, visited[1].second);
    }

    if (isWeightedDirectedGraph()) {
        ASSERT_EQ(2u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(1.0, visited[0].second);
        ASSERT_EQ(7.0, visited[1].second);
    }
}

TEST_P(GraphGTest, testForEdgesOf) {
    count m = 0;
    std::vector<int> visited(this->m_house, 0);

    this->Ghouse.forNodes([&](node u) {
        this->Ghouse.forEdgesOf(u, [&](node v, node w) {
            // edges should be v to w, so if we iterate over edges from u, u should be
            // equal v
            EXPECT_EQ(u, v);

            auto e = std::make_pair(v, w);
            auto it =
                std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
            if (!isDirected() && it == this->houseEdgesOut.end()) {
                auto e2 = std::make_pair(w, v);
                it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
                               e2);
            }

            EXPECT_TRUE(it != this->houseEdgesOut.end());

            // find index in edge array
            int i = std::distance(this->houseEdgesOut.begin(), it);
            if (isDirected()) {
                // make sure edge was not visited before (would be visited twice)
                EXPECT_EQ(0, visited[i]);
            }

            // mark edge as visited
            visited[i]++;
            m++;
        });
    });

    if (isDirected()) {
        // we iterated over all outgoing edges once
        EXPECT_EQ(this->m_house, m);
        for (auto c : visited) {
            EXPECT_EQ(1, c);
        }
    } else {
        // we iterated over all edges in both directions
        EXPECT_EQ(2 * this->m_house, m);
        for (auto c : visited) {
            EXPECT_EQ(2, c);
        }
    }
}

TEST_P(GraphGTest, testForWeightedEdgesOf) {
    count m = 0;
    std::vector<int> visited(this->m_house, 0);
    double sumOfWeights = 0;

    this->Ghouse.forNodes([&](node u) {
        this->Ghouse.forEdgesOf(u, [&](node v, node w, edgeweight ew) {
            // edges should be v to w, so if we iterate over edges from u, u should be
            // equal v
            EXPECT_EQ(u, v);
            sumOfWeights += ew;
            auto e = std::make_pair(v, w);
            auto it =
                std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
            if (!isDirected() && it == this->houseEdgesOut.end()) {
                auto e2 = std::make_pair(w, v);
                it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
                               e2);
            }

            EXPECT_TRUE(it != this->houseEdgesOut.end());

            // find index in edge array
            int i = std::distance(this->houseEdgesOut.begin(), it);
            if (isDirected()) {
                // make sure edge was not visited before (would be visited twice)
                EXPECT_EQ(0, visited[i]);
            }

            // mark edge as visited
            visited[i]++;
            m++;
        });
    });

    if (isGraph()) {
        EXPECT_EQ(sumOfWeights, m);
        EXPECT_EQ(2 * this->m_house, m);
        for (auto c : visited) {
            EXPECT_EQ(2, c);
        }
    }

    if (isWeightedGraph()) {
        // we iterated over all edges in both directions
        EXPECT_EQ(2 * this->m_house, m);
        EXPECT_EQ(sumOfWeights, 72);
        for (auto c : visited) {
            EXPECT_EQ(2, c);
        }
    }

    if (isDirectedGraph()) {
        // we iterated over all outgoing edges once
        EXPECT_EQ(this->m_house, m);
        EXPECT_EQ(sumOfWeights, m);
        for (auto c : visited) {
            EXPECT_EQ(1, c);
        }
    }

    if (isWeightedDirectedGraph()) {
        EXPECT_EQ(sumOfWeights, 36);
        EXPECT_EQ(this->m_house, m);
        for (auto c : visited) {
            EXPECT_EQ(1, c);
        }
    }
}

TEST_P(GraphGTest, testForInNeighborsOf) {
    std::vector<node> visited;
    this->Ghouse.forInNeighborsOf(2, [&](node v) { visited.push_back(v); });
    Aux::Parallel::sort(visited.begin(), visited.end());

    if (isDirected()) {
        EXPECT_EQ(2u, visited.size());
        EXPECT_EQ(0u, visited[0]);
        EXPECT_EQ(3u, visited[1]);
    } else {
        EXPECT_EQ(4u, visited.size());
        EXPECT_EQ(0u, visited[0]);
        EXPECT_EQ(1u, visited[1]);
        EXPECT_EQ(3u, visited[2]);
        EXPECT_EQ(4u, visited[3]);
    }
}

TEST_P(GraphGTest, testForWeightedInNeighborsOf) {
    std::vector<std::pair<node, edgeweight>> visited;
    this->Ghouse.forInNeighborsOf(3, [&](node v, edgeweight ew) {
        visited.push_back({v, ew});
    });
    Aux::Parallel::sort(visited.begin(), visited.end());

    if (isGraph()) {
        ASSERT_EQ(3u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(4u, visited[2].first);
        ASSERT_EQ(defaultEdgeWeight, visited[0].second);
        ASSERT_EQ(defaultEdgeWeight, visited[1].second);
        ASSERT_EQ(defaultEdgeWeight, visited[2].second);
    }

    if (isWeightedGraph()) {
        ASSERT_EQ(3u, visited.size());
        ASSERT_EQ(1u, visited[0].first);
        ASSERT_EQ(2u, visited[1].first);
        ASSERT_EQ(4u, visited[2].first);
        ASSERT_EQ(1.0, visited[0].second);
        ASSERT_EQ(7.0, visited[1].second);
        ASSERT_EQ(6.0, visited[2].second);
    }

    if (isDirectedGraph()) {
        ASSERT_EQ(1u, visited.size());
        ASSERT_EQ(4u, visited[0].first);
        ASSERT_EQ(defaultEdgeWeight, visited[0].second);
    }

    if (isWeightedDirectedGraph()) {
        ASSERT_EQ(1u, visited.size());
        ASSERT_EQ(4u, visited[0].first);
        ASSERT_EQ(6.0, visited[0].second);
    }
}

TEST_P(GraphGTest, forInEdgesOf) {
    std::vector<bool> visited(this->n_house, false);
    this->Ghouse.forInEdgesOf(3, [&](node u, node v) {
        ASSERT_EQ(3u, v);
        ASSERT_TRUE(this->Ahouse[u][v] > 0.0);
        ASSERT_TRUE(this->Ghouse.hasEdge(u, v));
        ASSERT_FALSE(visited[u]);
        visited[u] = true;
    });

    if (isDirected()) {
        EXPECT_FALSE(visited[0]);
        EXPECT_FALSE(visited[1]);
        EXPECT_FALSE(visited[2]);
        EXPECT_FALSE(visited[3]);
        EXPECT_TRUE(visited[4]);
    } else {
        EXPECT_FALSE(visited[0]);
        EXPECT_TRUE(visited[1]);
        EXPECT_TRUE(visited[2]);
        EXPECT_FALSE(visited[3]);
        EXPECT_TRUE(visited[4]);
    }
}

TEST_P(GraphGTest, testForWeightedInEdgesOf) {
    // add self-loop
    this->Ghouse.addEdge(3, 3, 2.5);
    this->Ahouse[3][3] = 2.5;

    std::vector<edgeweight> visited(this->n_house, -1.0);
    this->Ghouse.forInEdgesOf(3, [&](node v, node u, edgeweight ew) {
        ASSERT_EQ(-1.0, visited[v]);
        visited[u] = ew;
    });

    if (isGraph()) {
        ASSERT_EQ(-1.0, visited[0]);
        ASSERT_EQ(defaultEdgeWeight, visited[1]);
        ASSERT_EQ(defaultEdgeWeight, visited[2]);
        ASSERT_EQ(defaultEdgeWeight, visited[3]);
        ASSERT_EQ(defaultEdgeWeight, visited[4]);
    }

    if (isWeightedGraph()) {
        ASSERT_EQ(-1.0, visited[0]);
        ASSERT_EQ(this->Ahouse[3][1], visited[1]);
        ASSERT_EQ(this->Ahouse[3][2], visited[2]);
        ASSERT_EQ(this->Ahouse[3][3], visited[3]);
        ASSERT_EQ(this->Ahouse[3][4], visited[4]);
    }

    if (isDirectedGraph()) {
        ASSERT_EQ(-1.0, visited[0]);
        ASSERT_EQ(-1.0, visited[1]);
        ASSERT_EQ(-1.0, visited[2]);
        ASSERT_EQ(defaultEdgeWeight, visited[3]);
        ASSERT_EQ(defaultEdgeWeight, visited[4]);
    }

    if (isWeightedDirectedGraph()) {
        ASSERT_EQ(-1.0, visited[0]);
        ASSERT_EQ(-1.0, visited[1]);
        ASSERT_EQ(-1.0, visited[2]);
        ASSERT_EQ(this->Ahouse[3][3], visited[3]);
        ASSERT_EQ(this->Ahouse[4][3], visited[4]);
    }
}

/** REDUCTION ITERATORS **/

TEST_P(GraphGTest, testParallelSumForNodes) {
    count n = 10;
    Graph G = createGraph(n);
    double sum = G.parallelSumForNodes([](node v) { return 2 * v + 0.5; });

    double expected_sum = n * (n - 1) + n * 0.5;
    ASSERT_EQ(expected_sum, sum);
}

TEST_P(GraphGTest, testParallelSumForWeightedEdges) {
    double sum = this->Ghouse.parallelSumForEdges(
        [](node, node, edgeweight ew) { return 1.5 * ew; });

    double expected_sum = 1.5 * this->Ghouse.totalEdgeWeight();
    ASSERT_EQ(expected_sum, sum);
}

/** GRAPH SEARCHES **/

TEST_P(GraphGTest, testEdgeIndexGenerationDirected) {
    Graph G = Graph(10, false, true);
    G.addEdge(2, 0);
    G.addEdge(2, 1);
    G.addEdge(2, 2);
    G.addEdge(5, 6);
    G.addEdge(6, 5);
    G.addEdge(1, 2);

    G.indexEdges();

    // Check consecutiveness of edgeids according to edge iterators
    edgeid expectedId = 0;
    G.forEdges([&](node u, node v) { EXPECT_EQ(expectedId++, G.edgeId(u, v)); });

    // Add some more edges
    EXPECT_EQ(6, G.upperEdgeIdBound());
    G.addEdge(8, 9);
    EXPECT_EQ(7, G.upperEdgeIdBound());
    G.addEdge(9, 8);

    // Check that asymmetric edges do not have the same id
    EXPECT_NE(G.edgeId(6, 5), G.edgeId(5, 6));
    EXPECT_NE(G.edgeId(2, 1), G.edgeId(1, 2));
    EXPECT_NE(G.edgeId(9, 8), G.edgeId(8, 9));
    EXPECT_EQ(7, G.edgeId(9, 8));
    EXPECT_EQ(8, G.upperEdgeIdBound());
}

TEST_P(GraphGTest, testEdgeIndexGenerationUndirected) {
    Graph G = Graph(10, false, false);

    G.addEdge(0, 0);
    G.addEdge(2, 0);
    G.addEdge(2, 1);
    G.addEdge(2, 2);
    G.addEdge(5, 6);

    G.indexEdges();

    // Check consecutiveness of edgeids according to edge iterators
    edgeid expectedId = 0;
    G.forEdges([&](node u, node v) { EXPECT_EQ(expectedId++, G.edgeId(u, v)); });

    // Add some more edges. This will likely destroy consecutiveness...
    G.addEdge(3, 4);
    G.addEdge(7, 8);
    EXPECT_EQ(6, G.edgeId(7, 8));
    EXPECT_EQ(7, G.upperEdgeIdBound());

    // Anyway, check uniqueness and validity of the edgeids
    std::set<edgeid> ids;
    edgeid upperEdgeIdBound = G.upperEdgeIdBound();

    G.forEdges([&](node u, node v) {
        edgeid id = G.edgeId(u, v);
        EXPECT_EQ(id, G.edgeId(v, u));
        EXPECT_LT(id, upperEdgeIdBound);

        EXPECT_NE(none, id);
        EXPECT_FALSE(ids.erase(id));
        ids.insert(id);
    });
}

TEST_P(GraphGTest, testForEdgesWithIds) {
    std::vector<Graph> graphs;
    graphs.emplace_back(10, false, false);
    graphs.emplace_back(10, false, true);
    graphs.emplace_back(10, true, false);
    graphs.emplace_back(10, true, true);

    for (auto graph = graphs.begin(); graph != graphs.end(); ++graph) {
        graph->addEdge(0, 0);
        graph->addEdge(1, 2);
        graph->addEdge(4, 5);

        // No edge indices

        count m = 0;
        graph->forEdges([&](node, node, edgeid eid) {
            EXPECT_EQ(0, eid);
            m++;
        });
        ASSERT_EQ(3u, m);

        // With edge indices
        graph->indexEdges();

        edgeid expectedId = 0;
        m = 0;
        graph->forEdges([&](node, node, edgeid eid) {
            EXPECT_EQ(expectedId++, eid);
            EXPECT_LT(eid, graph->upperEdgeIdBound());
            m++;
        });
        ASSERT_EQ(3u, m);
    }
}

TEST_P(GraphGTest, testForWeightedEdgesWithIds) {
    std::vector<Graph> graphs;
    graphs.emplace_back(10, false, false);
    graphs.emplace_back(10, false, true);
    graphs.emplace_back(10, true, false);
    graphs.emplace_back(10, true, true);

    for (auto graph = graphs.begin(); graph != graphs.end(); ++graph) {
        graph->addEdge(0, 0, 2);
        graph->addEdge(1, 2, 2);
        graph->addEdge(4, 5, 2);

        // No edge indices

        count m = 0;
        edgeweight sum = 0;
        graph->forEdges([&](node, node, edgeweight ew, edgeid eid) {
            EXPECT_EQ(0, eid);
            m++;
            sum += ew;
        });
        ASSERT_EQ(3u, m);

        if (graph->isWeighted()) {
            ASSERT_EQ(6.0, sum);
        } else {
            ASSERT_EQ(3.0, sum);
        }

        // With edge indices
        graph->indexEdges();

        edgeid expectedId = 0;
        m = 0;
        sum = .0;
        graph->forEdges([&](node, node, edgeweight ew, edgeid eid) {
            EXPECT_EQ(expectedId++, eid);
            EXPECT_LT(eid, graph->upperEdgeIdBound());
            m++;
            sum += ew;
        });
        ASSERT_EQ(3u, m);

        if (graph->isWeighted()) {
            ASSERT_EQ(6.0, sum);
        } else {
            ASSERT_EQ(3.0, sum);
        }
    }
}

TEST_P(GraphGTest, testParallelForEdgesWithIds) {
    std::vector<Graph> graphs;
    graphs.emplace_back(10, false, false);
    graphs.emplace_back(10, false, true);
    graphs.emplace_back(10, true, false);
    graphs.emplace_back(10, true, true);

    for (auto graph = graphs.begin(); graph != graphs.end(); ++graph) {
        graph->addEdge(0, 0);
        graph->addEdge(1, 2);
        graph->addEdge(4, 5);

        // No edge indices
        count m = 0;
        edgeid sumedgeid = 0;
        graph->parallelForEdges([&](node, node, edgeid eid) {
#pragma omp atomic
            m++;
#pragma omp atomic
            sumedgeid += eid;
        });
        ASSERT_EQ(0, sumedgeid);
        ASSERT_EQ(3u, m);

        // With edge indices
        graph->indexEdges();

        m = 0;
        edgeid expectedId = 0;
        sumedgeid = 0;
        graph->parallelForEdges([&](node, node, edgeid eid) {
#pragma omp atomic
            expectedId++;
#pragma omp atomic
            sumedgeid += eid;
#pragma omp atomic
            m++;        
        });
        ASSERT_EQ(expectedId, graph->upperEdgeIdBound());
        ASSERT_EQ(sumedgeid, ((graph->upperEdgeIdBound()-1) * graph->upperEdgeIdBound()) / 2);
        ASSERT_EQ(3u, m);
    }
}

TEST_P(GraphGTest, testParallelForWeightedEdgesWithIds) {
    std::vector<Graph> graphs;
    graphs.emplace_back(10, false, false);
    graphs.emplace_back(10, false, true);
    graphs.emplace_back(10, true, false);
    graphs.emplace_back(10, true, true);

    for (auto graph = graphs.begin(); graph != graphs.end(); ++graph) {
        graph->addEdge(0, 0, 2);
        graph->addEdge(1, 2, 2);
        graph->addEdge(4, 5, 2);

        // No edge indices

        count m = 0;
        edgeweight sum = 0;
        edgeid sumedgeid = 0;
        graph->parallelForEdges([&](node, node, edgeweight ew, edgeid eid) {
#pragma omp atomic
            m++;
#pragma omp atomic
            sum += ew;
#pragma omp atomic
            sumedgeid += eid;
        });
        ASSERT_EQ(0, sumedgeid);
        ASSERT_EQ(3u, m);

        if (graph->isWeighted()) {
            ASSERT_EQ(6.0, sum);
        } else {
            ASSERT_EQ(3.0, sum);
        }

        // With edge indices
        graph->indexEdges();

        m = 0;
        sum = .0;
        edgeid expectedId = 0;
        sumedgeid = 0;
        graph->parallelForEdges([&](node, node, edgeweight ew, edgeid eid) {
#pragma omp atomic
            expectedId++;
#pragma omp atomic
            sumedgeid += eid;
#pragma omp atomic
            m++;
#pragma omp atomic
            sum += ew;
        });
        ASSERT_EQ(expectedId, graph->upperEdgeIdBound());
        ASSERT_EQ(sumedgeid, ((graph->upperEdgeIdBound()-1) * graph->upperEdgeIdBound()) / 2);
        ASSERT_EQ(3u, m);

        if (graph->isWeighted()) {
            ASSERT_EQ(6.0, sum);
        } else {
            ASSERT_EQ(3.0, sum);
        }
    }
}

/*TEST_P(GraphGTest, testInForEdgesUndirected) {
  METISGraphReader reader;
  Graph G = reader.read("input/PGPgiantcompo.graph");
  DEBUG(G.upperNodeIdBound());
  node u = 5474;
  G.forInEdgesOf(u, [&](node u, node z, edgeweight w){
    DEBUG("(1) node: ", u, " neigh:", z, " weight: ", w);
  });
  G.forEdgesOf(u, [&](node u, node z, edgeweight w){
    DEBUG("(2) node: ", u, " neigh:", z, " weight: ", w);
  });


  node source = 1492;
  DynBFS bfs(G, source, false);
  bfs.run();

  std::vector<std::pair<node, double> > choices1;
  G.forInEdgesOf(5474, [&](node t, node z, edgeweight w){
    INFO("considered edge (1): ", t, z, w);
    if (Aux::NumericTools::logically_equal(bfs.distance(t), bfs.distance(z) +
w)) {
      // workaround for integer overflow in large graphs
      bigfloat tmp = bfs.getNumberOfPaths(z) / bfs.getNumberOfPaths(t);
      double weight;
      tmp.ToDouble(weight);
      choices1.emplace_back(z, weight);
    }
  });
  std::vector<std::pair<node, double> > choices2;
  G.forEdgesOf(5474, [&](node t, node z, edgeweight w){
    INFO("considered edge (2): ", t, z, w);
    if (Aux::NumericTools::logically_equal(bfs.distance(t), bfs.distance(z) +
w)) {
      // workaround for integer overflow in large graphs
      bigfloat tmp = bfs.getNumberOfPaths(z) / bfs.getNumberOfPaths(t);
      double weight;
      tmp.ToDouble(weight);
      choices2.emplace_back(z, weight);
    }
  });

  INFO(choices1);
  INFO(choices2);
}
*/

TEST_P(GraphGTest, testCompactEdges) {
    Graph G = this->Ghouse;
    G.indexEdges();

    G.addEdge(0, 4);
    G.addEdge(0, 3);
    G.removeEdge(0, 3);
    G.removeEdge(0, 4);

    G.compactEdges();

    std::vector<std::pair<node, node>> outEdges;
    outEdges.reserve(this->Ghouse.numberOfEdges());

    this->Ghouse.forEdges([&](node u, node v) { outEdges.emplace_back(u, v); });

    auto it = outEdges.begin();

    G.forEdges([&](node u, node v) {
        ASSERT_NE(it, outEdges.end());
        EXPECT_EQ(it->first, u);
        EXPECT_EQ(it->second, v);
        ++it;
    });
}

TEST_P(GraphGTest, testSortEdges) {
    Graph G = this->Ghouse;

    for (int i = 0; i < 2; ++i) {
        if (i > 0) {
            G.indexEdges();
        }

        Graph origG = G;

        G.sortEdges();

        std::vector<std::tuple<node, node, edgeweight, edgeid>> edges;
        edges.reserve(origG.numberOfEdges() * 4);

        std::vector<std::tuple<node, edgeweight, edgeid>> outEdges;
        origG.forNodes([&](node u) {
            origG.forEdgesOf(u, [&](node, node v, edgeweight w, edgeid eid) {
                outEdges.emplace_back(v, w, eid);
            });

            Aux::Parallel::sort(outEdges.begin(), outEdges.end());

            for (auto x : outEdges) {
                edges.emplace_back(u, std::get<0>(x), std::get<1>(x), std::get<2>(x));
            }

            outEdges.clear();

            origG.forInEdgesOf(u, [&](node, node v, edgeweight w, edgeid eid) {
                outEdges.emplace_back(v, w, eid);
            });

            Aux::Parallel::sort(outEdges.begin(), outEdges.end());

            for (auto x : outEdges) {
                edges.emplace_back(u, std::get<0>(x), std::get<1>(x), std::get<2>(x));
            }

            outEdges.clear();
        });

        auto it = edges.begin();

        G.forNodes([&](node u) {
            G.forEdgesOf(u, [&](node u, node v, edgeweight w, edgeid eid) {
                ASSERT_NE(it, edges.end());
                EXPECT_EQ(*it, std::make_tuple(u, v, w, eid))
                    << "Out edge (" << u << ", " << v << ", " << w << ", " << eid
                    << ") was expected to be (" << std::get<0>(*it) << ", "
                    << std::get<1>(*it) << ", " << std::get<2>(*it) << ", "
                    << std::get<3>(*it) << ")";
                ++it;
            });
            G.forInEdgesOf(u, [&](node u, node v, edgeweight w, edgeid eid) {
                ASSERT_NE(it, edges.end());
                EXPECT_EQ(*it, std::make_tuple(u, v, w, eid))
                    << "In edge (" << u << ", " << v << ", " << w << ", " << eid
                    << ") was expected to be (" << std::get<0>(*it) << ", "
                    << std::get<1>(*it) << ", " << std::get<2>(*it) << ", "
                    << std::get<3>(*it) << ")";
                ++it;
            });
        });
    }
}

TEST_P(GraphGTest, testEdgeIdsAfterRemove) {
    constexpr node n = 100;

    Aux::Random::setSeed(42, true);
    auto G = createGraph(n, 10*n);
    G.indexEdges();
    auto original = G;

    // remove some nodes and edges
    G.removeNode(5);
    G.removeNode(10);
    while(2*G.numberOfEdges() > original.numberOfEdges()) {
        const auto e = GraphTools::randomEdge(G, false);
        G.removeEdge(e.first, e.second);
    }
    ASSERT_GT(G.numberOfEdges(), original.numberOfEdges() / 3);

    // check that the remaining edges still have the same ids
    G.forNodes([&] (node u) {
        G.forNeighborsOf(u, [&] (node, node v, edgeweight, edgeid id) {
            ASSERT_EQ(id, original.edgeId(u, v));
        });

        if (!isDirected()) return;

        G.forInNeighborsOf(u, [&] (node, node v, edgeweight, edgeid id) {
            ASSERT_EQ(id, original.edgeId(v, u));
        });
    });
}

} /* namespace NetworKit */
