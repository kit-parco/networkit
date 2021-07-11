// no-networkit-format
/*
 * GraphBuilderDirectSwapGTest.cpp
 *
 *  Created on: 14.08.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <tuple>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/Parallel.hpp>

#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/GraphBuilder.hpp>

namespace NetworKit {

class GraphBuilderDirectSwapGTest: public testing::TestWithParam< std::tuple<bool, bool> > {
public:
    virtual void SetUp();

protected:
    GraphBuilder bHouse;
    std::vector< std::pair<node, node> > houseEdgesOut;
    std::vector< std::vector<edgeweight> > Ahouse;
    count n_house;
    count m_house;

    bool isGraph() const { return !isWeighted() && !isDirected(); }
    bool isWeightedGraph() const { return isWeighted() && !isDirected(); }
    bool isDirectedGraph() const { return !isWeighted() && isDirected(); }
    bool isWeightedDirectedGraph() const { return isWeighted() && isDirected(); }

    bool isWeighted() const;
    bool isDirected() const;

    GraphBuilder createGraphBuilder(count n = 0) const;
    Graph toGraph(GraphBuilder& b) const;
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, GraphBuilderDirectSwapGTest, testing::Values(
    std::make_tuple(false, false),
    std::make_tuple(true, false),
    std::make_tuple(false, true),
    std::make_tuple(false, false)
));

bool GraphBuilderDirectSwapGTest::isWeighted() const {
    return std::get<0>(GetParam());
}
bool GraphBuilderDirectSwapGTest::isDirected() const {
    return std::get<1>(GetParam());
}

GraphBuilder GraphBuilderDirectSwapGTest::createGraphBuilder(count n) const {
    return GraphBuilder(n, isWeighted(), isDirected());
}

Graph GraphBuilderDirectSwapGTest::toGraph(GraphBuilder& b) const {
    return b.toGraph(false);
}

void GraphBuilderDirectSwapGTest::SetUp() {
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

    bHouse = createGraphBuilder(5);
    houseEdgesOut = {
        {3, 1},
        {1, 0},
        {0, 2},
        {2, 1},
        {1, 4},
        {4, 3},
        {3, 2},
        {2, 4}
    };
    Ahouse = {n_house, std::vector<edgeweight>(n_house, 0.0)};
    edgeweight ew = 1.0;
    for (auto& e : houseEdgesOut) {
        node u = e.first;
        node v = e.second;
        if (isDirected()) {
            bHouse.addHalfOutEdge(u, v, ew);
            bHouse.addHalfInEdge(v, u, ew);
        } else {
            bHouse.addHalfEdge(u, v, ew);
            bHouse.addHalfEdge(v, u, ew);
        }

        Ahouse[u][v] = ew;
        if (!bHouse.isDirected()) {
            Ahouse[v][u] = ew;
        }

        if (bHouse.isWeighted()) {
            ew += 1.0;
        }
    }
}

TEST_P(GraphBuilderDirectSwapGTest, testEmptyGraph) {
    auto b = createGraphBuilder();
    ASSERT_EQ(0u, b.numberOfNodes());

    Graph G = toGraph(b);

    ASSERT_EQ(0u, G.numberOfNodes());
    ASSERT_EQ(0u, G.numberOfEdges());
    ASSERT_TRUE(G.isEmpty());
}

TEST_P(GraphBuilderDirectSwapGTest, testAddNode) {
    auto b = createGraphBuilder();

    b.addNode();
    ASSERT_EQ(1u, b.numberOfNodes());

    b.addNode();
    b.addNode();
    ASSERT_EQ(3u, b.numberOfNodes());

    Graph G = toGraph(b);

    ASSERT_TRUE(G.hasNode(0));
    ASSERT_TRUE(G.hasNode(1));
    ASSERT_TRUE(G.hasNode(2));
    ASSERT_FALSE(G.hasNode(3));
    ASSERT_EQ(3u, G.numberOfNodes());
    ASSERT_EQ(0u, G.numberOfEdges());
    ASSERT_FALSE(G.isEmpty());
}


/** NODE PROPERTIES **/

TEST_P(GraphBuilderDirectSwapGTest, testDegree) {
    Graph Ghouse = toGraph(this->bHouse);
    if (isDirected()) {
        ASSERT_EQ(1u, Ghouse.degree(0));
        ASSERT_EQ(2u, Ghouse.degree(1));
        ASSERT_EQ(2u, Ghouse.degree(2));
        ASSERT_EQ(2u, Ghouse.degree(3));
        ASSERT_EQ(1u, Ghouse.degree(4));
    } else {
        ASSERT_EQ(2u, Ghouse.degree(0));
        ASSERT_EQ(4u, Ghouse.degree(1));
        ASSERT_EQ(4u, Ghouse.degree(2));
        ASSERT_EQ(3u, Ghouse.degree(3));
        ASSERT_EQ(3u, Ghouse.degree(4));
    }
}

TEST_P(GraphBuilderDirectSwapGTest, testDegreeIn) {
    Graph Ghouse = toGraph(this->bHouse);
    if (isDirected()) {
        ASSERT_EQ(1u, Ghouse.degreeIn(0));
        ASSERT_EQ(2u, Ghouse.degreeIn(1));
        ASSERT_EQ(2u, Ghouse.degreeIn(2));
        ASSERT_EQ(1u, Ghouse.degreeIn(3));
        ASSERT_EQ(2u, Ghouse.degreeIn(4));
    } else {
        ASSERT_EQ(2u, Ghouse.degreeIn(0));
        ASSERT_EQ(4u, Ghouse.degreeIn(1));
        ASSERT_EQ(4u, Ghouse.degreeIn(2));
        ASSERT_EQ(3u, Ghouse.degreeIn(3));
        ASSERT_EQ(3u, Ghouse.degreeIn(4));
    }
}

TEST_P(GraphBuilderDirectSwapGTest, testDegreeOut) {
    Graph Ghouse = toGraph(this->bHouse);
    if (isDirected()) {
        ASSERT_EQ(1u, Ghouse.degreeOut(0));
        ASSERT_EQ(2u, Ghouse.degreeOut(1));
        ASSERT_EQ(2u, Ghouse.degreeOut(2));
        ASSERT_EQ(2u, Ghouse.degreeOut(3));
        ASSERT_EQ(1u, Ghouse.degreeOut(4));
    } else {
        ASSERT_EQ(2u, Ghouse.degreeOut(0));
        ASSERT_EQ(4u, Ghouse.degreeOut(1));
        ASSERT_EQ(4u, Ghouse.degreeOut(2));
        ASSERT_EQ(3u, Ghouse.degreeOut(3));
        ASSERT_EQ(3u, Ghouse.degreeOut(4));
    }
}


/** EDGE MODIFIERS **/

TEST_P(GraphBuilderDirectSwapGTest, testAddHalfEdge) {
    auto b = createGraphBuilder(3);

    // Graph with 2 normal edges
    b.addHalfOutEdge(0, 1, 4.51);
    b.addHalfOutEdge(1, 2, 2.39);
    if (isDirected()) {
        b.addHalfInEdge(1, 0, 4.51);
        b.addHalfInEdge(2, 1, 2.39);
    } else {
        b.addHalfOutEdge(1, 0, 4.51);
        b.addHalfOutEdge(2, 1, 2.39);
    }

    Graph G = toGraph(b);

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
        // note: bidirectional edges are not supported, so both edges have different weights
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
}


/** GLOBAL PROPERTIES **/

TEST_P(GraphBuilderDirectSwapGTest, testIsWeighted) {
    ASSERT_EQ(isWeighted(), this->bHouse.isWeighted());
    Graph Ghouse = toGraph(this->bHouse);
    ASSERT_EQ(isWeighted(), Ghouse.isWeighted());
}

TEST_P(GraphBuilderDirectSwapGTest, testIsDirected) {
    ASSERT_EQ(isDirected(), this->bHouse.isDirected());
    Graph Ghouse = toGraph(this->bHouse);
    ASSERT_EQ(isDirected(), Ghouse.isDirected());
}

TEST_P(GraphBuilderDirectSwapGTest, testNumberOfSelfLoops) {
    auto b = createGraphBuilder(3);
    if (isDirected()) {
        b.addHalfOutEdge(0, 1);
        b.addHalfInEdge(1, 0);
        b.addHalfOutEdge(0, 0);
        b.addHalfInEdge(0, 0);
    } else {
        b.addHalfEdge(0, 1);
        b.addHalfEdge(0, 0);
        b.addHalfEdge(1,0);
    }
    Graph G = toGraph(b);
    ASSERT_EQ(1u, G.numberOfSelfLoops());
}

TEST_P(GraphBuilderDirectSwapGTest, testUpperNodeIdBound) {
    ASSERT_EQ(5u, this->bHouse.upperNodeIdBound());
    Graph Ghouse = toGraph(this->bHouse);
    ASSERT_EQ(5u, Ghouse.upperNodeIdBound());
}


/** EDGE ATTRIBUTES **/

TEST_P(GraphBuilderDirectSwapGTest, testSetWeight) {
    auto b = createGraphBuilder(10);
    b.addHalfOutEdge(0, 1);
    b.addHalfOutEdge(1, 2);
    if (isDirected()) {
        b.addHalfInEdge(1, 0);
        b.addHalfInEdge(2, 1);
    } else {
        b.addHalfOutEdge(1, 0);
        b.addHalfOutEdge(2, 1);
    }

    if (isWeighted()) {
        // edges should get weight defaultWeight on creation and setWeight should overwrite this
        b.setWeight(1, 2, 2.718);

        // setting an edge weight should create the edge if it doesn't exists

        b.setWeight(5, 6, 56.0);
        // directed graphs are not symmetric, undirected are
        b.setWeight(3, 4, 2.718);
        if (isDirected()) {
            b.setWeight(4, 3, 5.243);
        }

        // self-loop
        b.addHalfEdge(8, 8, 2.5);
        b.setWeight(8, 8, 3.14);

        b.setWeight(2, 1, 2.718);
        b.setWeight(6, 5, 56.0);
        b.setWeight(4, 3, 2.718);

        Graph G = toGraph(b);
        ASSERT_TRUE(G.checkConsistency());
        // edges should get weight defaultWeight on creation and setWeight should overwrite this
        ASSERT_EQ(defaultEdgeWeight, G.weight(0, 1));
        ASSERT_EQ(2.718, G.weight(1, 2));
        if (isDirected()) {
            ASSERT_EQ(nullWeight, G.weight(1, 0));
            ASSERT_EQ(nullWeight, G.weight(2, 1));
        } else {
            // undirected graph is symmetric
            ASSERT_EQ(defaultEdgeWeight, G.weight(1, 0));
            ASSERT_EQ(2.718, G.weight(2, 1));
        }

        // setting an edge weight should create the edge if it doesn't exists
        ASSERT_EQ(56.0, G.weight(5, 6));
        ASSERT_EQ(isDirected() ? nullWeight : 56.0, G.weight(6, 5));
        ASSERT_TRUE(G.hasEdge(5, 6));

        // directed graphs are not symmetric, undirected are
        if (isDirected()) {
            ASSERT_EQ(2.718, G.weight(3, 4));
            ASSERT_EQ(5.243, G.weight(4, 3));
        } else {
            // we have actually 2 edges in both direction. It is not defined which weight is returned.
            ASSERT_TRUE(G.weight(3, 4) == 2.718 || G.weight(3, 4) == 5.243);
            ASSERT_TRUE(G.weight(4, 3) == 2.718 || G.weight(4, 3) == 5.243);
        }

        // self-loop
        ASSERT_EQ(3.14, G.weight(8, 8));
    }
}


/** toGraph **/

TEST_P(GraphBuilderDirectSwapGTest, testSameAsGraph) {
    const double epsilon = 1e-6;
    const count runs = 10;
    const count n_max = 10;

    Aux::Random::setSeed(42, false);

    // in each loop run we will create a random graph (using a GraphBuilder and a Graph)
    // we will only use methods that both GraphBuilder and Graph support and
    for (index i = 0; i < runs; i++) {
        count n = Aux::Random::integer(n_max);
        auto b = createGraphBuilder(n);
        Graph G_expected(n, isWeighted(), isDirected());

        G_expected.forNodes([&](node v) {
            double p = Aux::Random::probability();
            // if we change edges we have to keep in mind, that GraphBuilder is designed to keep only half-edges.
            // e.g. if we have already added an edge v -> u, changing the weight of u -> v might create a new edge in the builder but change the existing edge in G_expected (does not apply for directed graphs)

            if (p < 0.1) { // new node
                n++;
                b.addNode();
                G_expected.addNode();
            } else { // new edge
                node u = Aux::Random::integer(v, n - 1); // self-loops possible
                edgeweight ew = Aux::Random::probability();
                G_expected.addEdge(v, u, ew);
                b.addHalfOutEdge(v, u, ew);
                if (isDirected()) {
                    b.addHalfInEdge(u, v, ew);
                } else if (u != v) {
                    b.addHalfOutEdge(u, v, ew);
                }
            }

            if (isWeighted()) {
                node u = Aux::Random::integer(v, n - 1); // self-loops possible
                edgeweight ew = Aux::Random::probability();
                if (p < 0.5) {
                    G_expected.setWeight(v, u, ew);
                    b.setOutWeight(v, u, ew);
                    if (isDirected()) {
                        b.setInWeight(u, v, ew);
                    } else if (u != v) {
                        b.setOutWeight(u, v, ew);
                    }
                } else {
                    G_expected.increaseWeight(v, u, ew);
                    b.increaseOutWeight(v, u, ew);
                    if (isDirected()) {
                        b.increaseInWeight(u, v, ew);
                    } else if (u != v) {
                        b.increaseOutWeight(u, v, ew);
                    }
                }
            }
        });

        Graph G_actual = toGraph(b);

        // check for correct graph properties
        ASSERT_EQ(G_expected.numberOfNodes(), G_actual.numberOfNodes());
        ASSERT_EQ(G_expected.numberOfEdges(), G_actual.numberOfEdges());
        ASSERT_EQ(G_expected.upperNodeIdBound(), G_actual.upperNodeIdBound());
        ASSERT_EQ(G_expected.numberOfSelfLoops(), G_actual.numberOfSelfLoops());

        // compare nodes and edges of G_expected and G_actual
        G_expected.forNodes([&](node v) {
            ASSERT_TRUE(G_actual.hasNode(v));
            ASSERT_EQ(G_expected.degree(v), G_actual.degree(v));
            ASSERT_EQ(G_expected.degreeIn(v), G_actual.degreeIn(v));
            ASSERT_EQ(G_expected.degreeOut(v), G_actual.degreeOut(v));
            ASSERT_NEAR(G_expected.weightedDegree(v), G_actual.weightedDegree(v), epsilon);
            ASSERT_NEAR(G_expected.weightedDegree(v, true),
                        G_actual.weightedDegree(v, true), epsilon);
        });
        G_expected.forEdges([&](node u, node v, edgeweight ew) {
            ASSERT_TRUE(G_actual.hasEdge(u, v));
            ASSERT_NEAR(ew, G_actual.weight(u, v), epsilon);
        });

        // make sure that G_actual has not more nodes/edges than G_expected
        G_actual.forNodes([&](node v) {
            ASSERT_TRUE(G_expected.hasNode(v));
        });
        G_actual.forEdges([&](node u, node v) {
            ASSERT_TRUE(G_expected.hasEdge(u, v));
        });
    }
}

TEST_P(GraphBuilderDirectSwapGTest, testForValidStateAfterToGraph) {
    Graph Ghouse = toGraph(this->bHouse);

    ASSERT_TRUE(this->bHouse.isEmpty());
    ASSERT_EQ(0u, this->bHouse.numberOfNodes());
    ASSERT_EQ(0u, this->bHouse.upperNodeIdBound());
    ASSERT_EQ(isWeighted(), this->bHouse.isWeighted());
    ASSERT_EQ(isDirected(), this->bHouse.isDirected());
    this->bHouse.forNodes([&](node) {
        FAIL();
    });

    Graph G1 = toGraph(this->bHouse);
    ASSERT_TRUE(G1.isEmpty());
    ASSERT_EQ(0u, G1.numberOfNodes());
    ASSERT_EQ(0u, G1.upperNodeIdBound());
    ASSERT_EQ(isWeighted(), G1.isWeighted());
    ASSERT_EQ(isDirected(), G1.isDirected());

    node v = this->bHouse.addNode();
    node u = this->bHouse.addNode();
    this->bHouse.addHalfOutEdge(v, u, 0.25);
    if (isDirected()) {
        this->bHouse.addHalfInEdge(u, v, 0.25);
    } else {
        this->bHouse.addHalfOutEdge(u, v, 0.25);
    }

    Graph G2 = toGraph(this->bHouse);
    ASSERT_FALSE(G2.isEmpty());
    ASSERT_EQ(2u, G2.numberOfNodes());
    ASSERT_EQ(1u, G2.numberOfEdges());
    ASSERT_TRUE(G2.hasEdge(v, u));
    if (!isDirected()) {
        ASSERT_TRUE(G2.hasEdge(u, v));
    }
    if (isWeighted()) {
        ASSERT_EQ(0.25, G2.weight(v, u));
    } else {
        ASSERT_EQ(1.0, G2.weight(v, u));
    }
}

/** NODE ITERATORS **/

TEST_P(GraphBuilderDirectSwapGTest, testForNodes) {
    auto b = createGraphBuilder(3);
    std::vector<bool> visited(4, false);
    b.forNodes([&](node v) {
        ASSERT_FALSE(visited[v]);
        if (v == 2) {
            b.addNode();
        }
        visited[v] = true;
    });
    for (bool x : visited) {
        ASSERT_TRUE(x);
    }
}

TEST_P(GraphBuilderDirectSwapGTest, testParallelForNodes) {
    std::vector<node> visited(bHouse.upperNodeIdBound());
    this->bHouse.parallelForNodes([&](node u) {
        visited[u] = u;
    });

    Aux::Parallel::sort(visited.begin(), visited.end());

    ASSERT_EQ(5u, visited.size());
    for (index i = 0; i < this->bHouse.upperNodeIdBound(); i++) {
        ASSERT_EQ(i, visited[i]);
    }
}

TEST_P(GraphBuilderDirectSwapGTest, testForNodePairs) {
    count n = 10;
    auto b = createGraphBuilder(n);

    std::vector< std::vector<bool> > visited(n, std::vector<bool>(n, false));
    b.forNodePairs([&](node u, node v) {
        if (visited[u][v] || visited[v][u]) {
            FAIL();
        } else {
            visited[u][v] = true;
        }
    });

    for (node u = 0; u < n; u++) {
        for (node v = 0; v < n; v++) {
            if (u == v) {
                ASSERT_FALSE(visited[u][u]);
            } else {
                ASSERT_TRUE(visited[u][v] ^ visited[v][u]);
            }
        }
    }
}

TEST_P(GraphBuilderDirectSwapGTest, testParallelForNodePairs) {
    count n = 10;
    auto b = createGraphBuilder(n);

    std::vector< std::vector<bool> > visited(n, std::vector<bool>(n, false));
    b.forNodePairs([&](node u, node v) {
        if (visited[u][v]) {
            FAIL();
        } else {
            visited[u][v] = true;
        }
    });

    for (node u = 0; u < n; u++) {
        for (node v = 0; v < n; v++) {
            if (u == v) {
                ASSERT_FALSE(visited[u][u]);
            } else {
                ASSERT_TRUE(visited[u][v] ^ visited[v][u]);
            }
        }
    }
}

} /* namespace NetworKit */
