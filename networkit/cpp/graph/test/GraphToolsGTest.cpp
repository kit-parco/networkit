/*
 * GraphToolsGTest.cpp
 *
 *  Created on: 22.11.14
 *      Author: Maximilian Vogel
 */

#include <array>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class GraphToolsGTest : public testing::TestWithParam<std::pair<bool, bool>> {
protected:
    Graph generateRandomWeights(const Graph &G) const;
    bool weighted() const noexcept;
    bool directed() const noexcept;
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, GraphToolsGTest,
                         testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                                         std::make_pair(false, true), std::make_pair(true, true)));

Graph GraphToolsGTest::generateRandomWeights(const Graph &G) const {
    Graph Gw(G, true, G.isDirected());
    Gw.forEdges([&](node u, node v) { Gw.setWeight(u, v, Aux::Random::probability()); });
    return Gw;
}

bool GraphToolsGTest::weighted() const noexcept {
    return GetParam().first;
}

bool GraphToolsGTest::directed() const noexcept {
    return GetParam().second;
}

TEST_P(GraphToolsGTest, testSize) {
    constexpr count n = 100;
    constexpr double p = 0.1;
    constexpr count updates = 10;

    auto doTest = [](const Graph &G) {
        const auto size = GraphTools::size(G);
        EXPECT_EQ(size.first, G.numberOfNodes());
        EXPECT_EQ(size.second, G.numberOfEdges());
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(n, p, directed()).generate();
        if (weighted()) {
            G = generateRandomWeights(G);
        }

        doTest(G);
        for (count i = 0; i < updates; ++i) {
            G.removeNode(GraphTools::randomNode(G));
            doTest(G);
        }

        for (count i = 0; i < updates && G.numberOfEdges(); ++i) {
            const auto randomEdge = GraphTools::randomEdge(G);
            G.removeEdge(randomEdge.first, randomEdge.second);
            doTest(G);
        }
    }
}

TEST_P(GraphToolsGTest, testDensity) {
    constexpr count n = 100;
    constexpr double p = 0.1;

    auto doTest = [](const Graph &G) {
        const auto density = GraphTools::density(G);
        EXPECT_TRUE(density >= 0);
        EXPECT_EQ(density > 0, G.numberOfNodes() > 1 && G.numberOfEdges() - G.numberOfSelfLoops());
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(n, p, directed()).generate();
        if (weighted()) {
            G = generateRandomWeights(G);
        }

        doTest(G);
        for (node u = 0; u < G.upperNodeIdBound(); ++u) {
            if (G.hasNode(u))
                G.removeNode(u);
            doTest(G);
        }
    }
}

TEST_P(GraphToolsGTest, testVolume) {
    constexpr count n = 100;
    constexpr double p = 0.1;

    auto doTest = [&](const Graph &G) {
        // Volume for graph G
        const auto volume = GraphTools::volume(G);

        // Volume for either directed/undirected
        double mod = directed() ? 1.0 : 2.0;
        double baseG = weighted() ? G.totalEdgeWeight() : static_cast<double>(G.numberOfEdges());
        EXPECT_DOUBLE_EQ(volume, mod * baseG);

        // Volume for subgraph G2
        node seed = Aux::Random::integer(G.upperNodeIdBound());

        const auto first = G.nodeRange().begin();
        const auto lastOfSet = G.nodeRange().begin().operator++(seed);
        const auto last = G.nodeRange().end();

        const auto volumeSet = GraphTools::volume(G, first, lastOfSet);
        const auto inVolumeSet = GraphTools::inVolume(G, first, lastOfSet);
        const auto volumeComplementSet = GraphTools::volume(G, lastOfSet, last);
        const auto inVolumeComplementSet = GraphTools::inVolume(G, lastOfSet, last);

        EXPECT_NEAR(volumeSet + volumeComplementSet, volume, 1e-7);
        EXPECT_NEAR(inVolumeSet + inVolumeComplementSet, volume, 1e-7);
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(n, p, directed()).generate();
        if (weighted()) {
            G = generateRandomWeights(G);
        }

        doTest(G);
    }
}

TEST_P(GraphToolsGTest, testMaxDegree) {
    constexpr count n = 100;
    constexpr double p = 0.1;
    constexpr count edgeUpdates = 10;

    auto computeMaxDeg = [&](const Graph &G, bool inDegree) {
        count maxDeg = 0;
        G.forNodes([&](const node u) {
            maxDeg = std::max(maxDeg, inDegree ? G.degreeIn(u) : G.degreeOut(u));
        });

        return maxDeg;
    };

    auto computeMaxWeightedDeg = [&](const Graph &G, bool inDegree) {
        edgeweight maxDeg = std::numeric_limits<edgeweight>::min();
        G.forNodes([&](const node u) {
            maxDeg = std::max(maxDeg, inDegree ? G.weightedDegreeIn(u) : G.weightedDegree(u));
        });

        return maxDeg;
    };

    auto doTest = [&](const Graph &G) {
        EXPECT_EQ(GraphTools::maxDegree(G), computeMaxDeg(G, false));
        EXPECT_EQ(GraphTools::maxInDegree(G), computeMaxDeg(G, true));
        EXPECT_DOUBLE_EQ(GraphTools::maxWeightedDegree(G), computeMaxWeightedDeg(G, false));
        EXPECT_DOUBLE_EQ(GraphTools::maxWeightedInDegree(G), computeMaxWeightedDeg(G, true));
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(n, p, directed()).generate();
        if (weighted()) {
            G = generateRandomWeights(G);
        }

        doTest(G);
        for (count i = 0; i < edgeUpdates; ++i) {
            const auto e = GraphTools::randomEdge(G);
            G.removeEdge(e.first, e.second);
            doTest(G);
        }

        for (count i = 0; i < edgeUpdates; ++i) {
            node u = GraphTools::randomNode(G);
            node v = GraphTools::randomNode(G);
            while (G.hasEdge(u, v)) {
                u = GraphTools::randomNode(G);
                v = GraphTools::randomNode(G);
            }
            G.addEdge(u, v);
            doTest(G);
        }
    }
}

TEST_P(GraphToolsGTest, testRandomNode) {
    constexpr count n = 20;
    constexpr double p = 0.01;
    constexpr count samples = 100000;
    constexpr double maxAbsoluteError = 0.005;

    Aux::Random::setSeed(42, false);

    auto G = ErdosRenyiGenerator(n, p, directed()).generate();
    if (weighted()) {
        G = Graph(G, true, G.isDirected());
    }

    std::vector<count> drawCounts(n, 0);
    for (count i = 0; i < samples; i++) {
        ++drawCounts[GraphTools::randomNode(G)];
    }

    for (node v = 0; v < n; v++) {
        const auto p = static_cast<double>(drawCounts[v]) / static_cast<double>(samples);
        ASSERT_NEAR(1.0 / n, p, maxAbsoluteError);
    }
}

TEST_P(GraphToolsGTest, testRandomNodes) {
    const count n = 20;
    const double p = 0.01;

    Aux::Random::setSeed(42, false);

    auto G = ErdosRenyiGenerator(n, p, directed()).generate();
    if (weighted()) {
        G = Graph(G, true, G.isDirected());
    }

    count sampleAll = G.numberOfNodes();
    std::vector<node> allSample = GraphTools::randomNodes(G, sampleAll);
    EXPECT_EQ(allSample.size(), sampleAll);
    EXPECT_EQ((std::unordered_set<node>{allSample.begin(), allSample.end()}).size(), sampleAll);

    count sampleMost = 15;
    std::vector<node> mostSample = GraphTools::randomNodes(G, sampleMost);
    EXPECT_EQ(mostSample.size(), sampleMost);
    EXPECT_EQ((std::unordered_set<node>{mostSample.begin(), mostSample.end()}).size(), sampleMost);

    count sampleSome = 5;
    std::vector<node> someSample = GraphTools::randomNodes(G, sampleSome);
    EXPECT_EQ(someSample.size(), sampleSome);
    EXPECT_EQ((std::unordered_set<node>{someSample.begin(), someSample.end()}).size(), sampleSome);
}

TEST_P(GraphToolsGTest, testRandomNeighbor) {
    constexpr count n = 10;
    auto G = Graph(n, weighted(), directed());
    G.addEdge(2, 0);
    G.addEdge(2, 1);
    G.addEdge(2, 2);
    G.addEdge(5, 6);

    Aux::Random::setSeed(42, false);

    ASSERT_EQ(none, GraphTools::randomNeighbor(G, 3));
    ASSERT_EQ(6u, GraphTools::randomNeighbor(G, 5));

    if (G.isDirected()) {
        ASSERT_EQ(none, GraphTools::randomNeighbor(G, 1));
    } else {
        ASSERT_EQ(2u, GraphTools::randomNeighbor(G, 1));
    }

    constexpr count nn = 3;
    constexpr count samples = 100000;
    constexpr double maxAbsoluteError = 0.005;
    std::vector<count> drawCounts(nn, 0);
    for (count i = 0; i < samples; i++) {
        ++drawCounts[GraphTools::randomNeighbor(G, 2)];
    }
    for (node v = 0; v < nn; v++) {
        double p = static_cast<double>(drawCounts[v]) / static_cast<double>(samples);
        ASSERT_NEAR(1.0 / nn, p, maxAbsoluteError);
    }
}

TEST_P(GraphToolsGTest, testRandomEdge) {
    Aux::Random::setSeed(1, false);
    // we only test the uniform version
    constexpr count n = 4;
    constexpr count m = 5;
    constexpr count samples = 100000;
    constexpr double maxAbsoluteError = 0.005;

    Graph G(n, weighted(), directed());
    G.addEdge(0, 1); // 0 * 1 = 0
    G.addEdge(1, 2); // 1 * 2 = 2
    G.addEdge(3, 2); // 3 * 2 = 1 (mod 5)
    G.addEdge(2, 2); // 2 * 2 = 4
    G.addEdge(3, 1); // 3 * 1 = 3
    ASSERT_EQ(m, G.numberOfEdges());

    std::vector<count> drawCounts(m, 0);
    for (count i = 0; i < samples; ++i) {
        const auto e = GraphTools::randomEdge(G, true);
        count id = (e.first * e.second) % 5;
        drawCounts[id]++;
    }
    for (node id = 0; id < m; id++) {
        double p = static_cast<double>(drawCounts[id]) / static_cast<double>(samples);
        ASSERT_NEAR(1.0 / m, p, maxAbsoluteError);
    }
}

TEST_P(GraphToolsGTest, testRandomEdges) {
    Aux::Random::setSeed(1, false);
    // we only test the uniform version
    constexpr count n = 4;
    constexpr count m = 5;
    constexpr count samples = 100000;
    constexpr double maxAbsoluteError = 0.005;

    Graph G(n, weighted(), directed());
    G.addEdge(0, 1); // 0 * 1 = 0
    G.addEdge(1, 2); // 1 * 2 = 2
    G.addEdge(3, 2); // 3 * 2 = 1 (mod 5)
    G.addEdge(2, 2); // 2 * 2 = 4
    G.addEdge(3, 1); // 3 * 1 = 3
    ASSERT_EQ(m, G.numberOfEdges());

    std::vector<count> drawCounts(m, 0);
    for (const auto &e : GraphTools::randomEdges(G, samples)) {
        const auto id = (e.first * e.second) % 5;
        ++drawCounts[id];
    }
    for (node id = 0; id < m; id++) {
        const double p = static_cast<double>(drawCounts[id]) / static_cast<double>(samples);
        ASSERT_NEAR(1.0 / static_cast<double>(m), p, maxAbsoluteError);
    }
}

TEST_P(GraphToolsGTest, testGetContinuousOnContinuous) {
    Graph G(10, weighted(), directed());
    auto nodeIds = GraphTools::getContinuousNodeIds(G);
    std::unordered_map<node, node> reference = {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4},
                                                {5, 5}, {6, 6}, {7, 7}, {8, 8}, {9, 9}};
    EXPECT_EQ(reference, nodeIds);
}

TEST_P(GraphToolsGTest, testGetContinuousOnDeletedNodes1) {
    Graph G(10, weighted(), directed());
    G.removeNode(0);
    G.removeNode(1);
    G.removeNode(2);
    G.removeNode(3);
    G.removeNode(4);
    auto nodeIds = GraphTools::getContinuousNodeIds(G);
    std::unordered_map<node, node> reference = {{5, 0}, {6, 1}, {7, 2}, {8, 3}, {9, 4}};
    EXPECT_EQ(reference, nodeIds);
}

TEST_P(GraphToolsGTest, testGetContinuousOnDeletedNodes2) {
    Graph G(10, weighted(), directed());
    G.removeNode(0);
    G.removeNode(2);
    G.removeNode(4);
    G.removeNode(6);
    G.removeNode(8);
    auto nodeIds = GraphTools::getContinuousNodeIds(G);
    std::unordered_map<node, node> reference = {{1, 0}, {3, 1}, {5, 2}, {7, 3}, {9, 4}};
    EXPECT_EQ(reference, nodeIds);
}

TEST_F(GraphToolsGTest, testGetCompactedGraphUndirectedUnweighted1) {
    Graph G(10, false, false);
    G.addEdge(0, 1);
    G.addEdge(2, 1);
    G.addEdge(0, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 6);
    G.addEdge(4, 8);
    G.addEdge(5, 9);
    G.addEdge(3, 7);
    G.addEdge(5, 7);

    auto nodeMap = GraphTools::getContinuousNodeIds(G);
    auto Gcompact = GraphTools::getCompactedGraph(G, nodeMap);

    EXPECT_EQ(G.numberOfNodes(), Gcompact.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), Gcompact.numberOfEdges());
    EXPECT_EQ(G.isDirected(), Gcompact.isDirected());
    EXPECT_EQ(G.isWeighted(), Gcompact.isWeighted());
    // TODOish: find a deeper test to check if the structure of the graphs are the same,
    // probably compare results of some algorithms or compare each edge with a reference node id
    // map.
}

TEST_F(GraphToolsGTest, testGetCompactedGraphUndirectedUnweighted2) {
    Graph G(10, false, false);
    G.removeNode(0);
    G.removeNode(2);
    G.removeNode(4);
    G.removeNode(6);
    G.removeNode(8);
    G.addEdge(1, 3);
    G.addEdge(5, 3);
    G.addEdge(7, 5);
    G.addEdge(7, 9);
    G.addEdge(1, 9);

    auto nodeMap = GraphTools::getContinuousNodeIds(G);
    auto Gcompact = GraphTools::getCompactedGraph(G, nodeMap);

    EXPECT_NE(G.upperNodeIdBound(), Gcompact.upperNodeIdBound());
    EXPECT_EQ(G.numberOfNodes(), Gcompact.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), Gcompact.numberOfEdges());
    EXPECT_EQ(G.isDirected(), Gcompact.isDirected());
    EXPECT_EQ(G.isWeighted(), Gcompact.isWeighted());
    // TODOish: find a deeper test to check if the structure of the graphs are the same,
    // probably compare results of some algorithms or compare each edge with a reference node id
    // map.
}

TEST_F(GraphToolsGTest, testGetCompactedGraphUndirectedWeighted1) {
    Graph G(10, true, false);
    G.removeNode(0);
    G.removeNode(2);
    G.removeNode(4);
    G.removeNode(6);
    G.removeNode(8);
    G.addEdge(1, 3, 0.2);
    G.addEdge(5, 3, 2132.351);
    G.addEdge(7, 5, 3.14);
    G.addEdge(7, 9, 2.7);
    G.addEdge(1, 9, 0.12345);

    auto nodeMap = GraphTools::getContinuousNodeIds(G);
    auto Gcompact = GraphTools::getCompactedGraph(G, nodeMap);

    EXPECT_EQ(G.totalEdgeWeight(), Gcompact.totalEdgeWeight());
    EXPECT_NE(G.upperNodeIdBound(), Gcompact.upperNodeIdBound());
    EXPECT_EQ(G.numberOfNodes(), Gcompact.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), Gcompact.numberOfEdges());
    EXPECT_EQ(G.isDirected(), Gcompact.isDirected());
    EXPECT_EQ(G.isWeighted(), Gcompact.isWeighted());
    // TODOish: find a deeper test to check if the structure of the graphs are the same,
    // probably compare results of some algorithms or compare each edge with a reference node id
    // map.
}

TEST_F(GraphToolsGTest, testGetCompactedGraphDirectedWeighted1) {
    Graph G(10, true, true);
    G.removeNode(0);
    G.removeNode(2);
    G.removeNode(4);
    G.removeNode(6);
    G.removeNode(8);
    G.addEdge(1, 3, 0.2);
    G.addEdge(5, 3, 2132.351);
    G.addEdge(7, 5, 3.14);
    G.addEdge(7, 9, 2.7);
    G.addEdge(1, 9, 0.12345);

    auto nodeMap = GraphTools::getContinuousNodeIds(G);
    auto Gcompact = GraphTools::getCompactedGraph(G, nodeMap);

    EXPECT_EQ(G.totalEdgeWeight(), Gcompact.totalEdgeWeight());
    EXPECT_NE(G.upperNodeIdBound(), Gcompact.upperNodeIdBound());
    EXPECT_EQ(G.numberOfNodes(), Gcompact.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), Gcompact.numberOfEdges());
    EXPECT_EQ(G.isDirected(), Gcompact.isDirected());
    EXPECT_EQ(G.isWeighted(), Gcompact.isWeighted());
    // TODOish: find a deeper test to check if the structure of the graphs are the same,
    // probably compare results of some algorithms or compare each edge with a reference node id
    // map.
}

TEST_F(GraphToolsGTest, testGetCompactedGraphDirectedUnweighted1) {
    Graph G(10, false, true);
    G.removeNode(0);
    G.removeNode(2);
    G.removeNode(4);
    G.removeNode(6);
    G.removeNode(8);
    G.addEdge(1, 3);
    G.addEdge(5, 3);
    G.addEdge(7, 5);
    G.addEdge(7, 9);
    G.addEdge(1, 9);
    auto nodeMap = GraphTools::getContinuousNodeIds(G);
    auto Gcompact = GraphTools::getCompactedGraph(G, nodeMap);

    EXPECT_EQ(G.totalEdgeWeight(), Gcompact.totalEdgeWeight());
    EXPECT_NE(G.upperNodeIdBound(), Gcompact.upperNodeIdBound());
    EXPECT_EQ(G.numberOfNodes(), Gcompact.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), Gcompact.numberOfEdges());
    EXPECT_EQ(G.isDirected(), Gcompact.isDirected());
    EXPECT_EQ(G.isWeighted(), Gcompact.isWeighted());
    // TODOish: find a deeper test to check if the structure of the graphs are the same,
    // probably compare results of some algorithms or compare each edge with a reference node id
    // map.
}

TEST_P(GraphToolsGTest, testInvertedMapping) {
    Graph G(10, weighted(), directed());
    G.removeNode(0);
    G.removeNode(2);
    G.removeNode(4);
    G.removeNode(6);
    G.removeNode(8);
    G.addEdge(1, 3);
    G.addEdge(5, 3);
    G.addEdge(7, 5);
    G.addEdge(7, 9);
    G.addEdge(1, 9);
    auto nodeMap = GraphTools::getContinuousNodeIds(G);
    auto invertedNodeMap = GraphTools::invertContinuousNodeIds(nodeMap, G);

    EXPECT_EQ(6, invertedNodeMap.size());

    std::vector<node> reference = {1, 3, 5, 7, 9, 10};
    EXPECT_EQ(reference, invertedNodeMap);
}

TEST_F(GraphToolsGTest, testRestoreGraph) {
    Graph G(10, false, true);
    G.removeNode(0);
    G.removeNode(2);
    G.removeNode(4);
    G.removeNode(6);
    G.removeNode(8);
    G.addEdge(1, 3);
    G.addEdge(5, 3);
    G.addEdge(7, 5);
    G.addEdge(7, 9);
    G.addEdge(1, 9);
    auto nodeMap = GraphTools::getContinuousNodeIds(G);
    auto invertedNodeMap = GraphTools::invertContinuousNodeIds(nodeMap, G);
    std::vector<node> reference = {1, 3, 5, 7, 9, 10};

    EXPECT_EQ(6, invertedNodeMap.size());
    EXPECT_EQ(reference, invertedNodeMap);

    auto Gcompact = GraphTools::getCompactedGraph(G, nodeMap);
    Graph Goriginal = GraphTools::restoreGraph(invertedNodeMap, Gcompact);

    EXPECT_EQ(Goriginal.totalEdgeWeight(), Gcompact.totalEdgeWeight());
    EXPECT_NE(Goriginal.upperNodeIdBound(), Gcompact.upperNodeIdBound());
    EXPECT_EQ(Goriginal.numberOfNodes(), Gcompact.numberOfNodes());
    EXPECT_EQ(Goriginal.numberOfEdges(), Gcompact.numberOfEdges());
    EXPECT_EQ(Goriginal.isDirected(), Gcompact.isDirected());
    EXPECT_EQ(Goriginal.isWeighted(), Gcompact.isWeighted());
}

TEST_P(GraphToolsGTest, testGetRemappedGraph) {
    const auto n = 4;
    Graph G(n, weighted(), directed());
    for (auto i : {0, 1, 2})
        G.addEdge(i, i + 1, i);

    if (directed())
        G.addEdge(1, 1, 12);

    std::vector<node> perm(n);
    for (int i = 0; i < n; ++i)
        perm[i] = i;

    std::mt19937_64 gen;
    for (int iter = 0; iter < 10; iter++) {
        std::shuffle(perm.begin(), perm.end(), gen);
        auto G1 = GraphTools::getRemappedGraph(G, n, [&](node i) { return perm[i]; });
        ASSERT_EQ(G1.numberOfNodes(), n);
        ASSERT_EQ(G1.numberOfEdges(), G.numberOfEdges());
        ASSERT_EQ(G1.numberOfSelfLoops(), G.numberOfSelfLoops());

        for (int i = 0; i < n; ++i) {
            for (int j = 0; i < n; ++i) {
                ASSERT_EQ(G.hasEdge(i, j), G1.hasEdge(perm[i], perm[j]));
                ASSERT_EQ(G.weight(i, j), G1.weight(perm[i], perm[j]));
            }
        }
    }
}

TEST_P(GraphToolsGTest, testGetRemappedGraphWithDelete) {
    const auto n = 4;
    Graph G(n, weighted(), directed());
    for (auto i : {0, 1, 2})
        G.addEdge(i, i + 1, i);

    if (directed())
        G.addEdge(1, 1, 12);

    std::vector<node> perm(n);
    for (int i = 0; i < n; ++i)
        perm[i] = i;

    std::mt19937_64 gen;
    std::uniform_int_distribution<node> distr(0, n - 1);
    for (int iter = 0; iter < 10; iter++) {
        std::shuffle(perm.begin(), perm.end(), gen);

        const auto del = distr(gen);

        auto G1 = GraphTools::getRemappedGraph(
            G, n, [&](node i) { return perm[i]; }, [&](node i) { return i == del; });

        auto expected_num_edges = G.numberOfEdges();
        expected_num_edges -= G.degree(del);
        if (directed())
            expected_num_edges -= G.degreeIn(del);
        // do double count self-loops
        expected_num_edges += G.hasEdge(del, del);

        ASSERT_EQ(G1.numberOfNodes(), n);
        ASSERT_EQ(G1.numberOfEdges(), expected_num_edges) << " del=" << del;
        ASSERT_EQ(G1.numberOfSelfLoops(), G.numberOfSelfLoops() - G.hasEdge(del, del))
            << " del=" << del;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; i < n; ++i) {
                if (i == static_cast<int>(del) || j == static_cast<int>(del)) {
                    ASSERT_FALSE(G1.hasEdge(perm[i], perm[j]))
                        << "i=" << i << " j=" << j << " del=" << del;
                } else {
                    ASSERT_EQ(G.hasEdge(i, j), G1.hasEdge(perm[i], perm[j]));
                    ASSERT_EQ(G.weight(i, j), G1.weight(perm[i], perm[j]));
                }
            }
        }
    }
}

TEST_P(GraphToolsGTest, testCopyNodes) {
    constexpr count n = 200;
    constexpr double p = 0.01;
    constexpr count nodesToDelete = 50;

    auto checkNodes = [&](const Graph &G, const Graph &GCopy) {
        EXPECT_EQ(G.isDirected(), GCopy.isDirected());
        EXPECT_EQ(G.isWeighted(), GCopy.isWeighted());
        EXPECT_EQ(G.numberOfNodes(), GCopy.numberOfNodes());
        EXPECT_EQ(GCopy.numberOfEdges(), 0);
        for (node u = 0; u < G.upperNodeIdBound(); ++u) {
            EXPECT_EQ(G.hasNode(u), GCopy.hasNode(u));
        }
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(n, p, directed()).generate();

        auto GCopy = GraphTools::copyNodes(G);
        checkNodes(G, GCopy);
        for (count i = 0; i < nodesToDelete; ++i) {
            G.removeNode(GraphTools::randomNode(G));
            GCopy = GraphTools::copyNodes(G);
            checkNodes(G, GCopy);
        }
    }
}

TEST_P(GraphToolsGTest, testSubgraphFromNodesUndirected) {
    auto G = Graph(4, weighted(), false);

    /**
     *      1
     *   /  |  \
     * 0    |    3
     *   \  |  /
     *      2
     */

    G.addEdge(0, 1, 1.0);
    G.addEdge(0, 2, 2.0);
    G.addEdge(3, 1, 4.0);
    G.addEdge(3, 2, 5.0);
    G.addEdge(1, 2, 3.0);

    std::array<bool, 2> compactOptions{{true, false}};

    for (bool compact : compactOptions) {
        std::unordered_set<node> nodes = {0};
        auto res = GraphTools::subgraphFromNodes(G, nodes.begin(), nodes.end(), compact);
        EXPECT_EQ(weighted(), res.isWeighted());
        EXPECT_FALSE(res.isDirected());
        EXPECT_EQ(res.numberOfNodes(), 1);
        EXPECT_EQ(res.numberOfEdges(), 0);
    }

    {
        std::unordered_set<node> nodes = {0};
        auto res = GraphTools::subgraphAndNeighborsFromNodes(G, nodes, true);

        EXPECT_EQ(res.numberOfNodes(), 3);
        EXPECT_EQ(res.numberOfEdges(), 2); // 0-1, 0-2, NOT 1-2

        EXPECT_DOUBLE_EQ(G.weight(0, 1), weighted() ? 1.0 : defaultEdgeWeight);
        EXPECT_DOUBLE_EQ(G.weight(0, 2), weighted() ? 2.0 : defaultEdgeWeight);
    }

    for (bool compact : compactOptions) {
        std::unordered_set<node> nodes = {0, 1};
        auto res = GraphTools::subgraphFromNodes(G, nodes.begin(), nodes.end(), compact);
        EXPECT_EQ(res.numberOfNodes(), 2);
        EXPECT_EQ(res.numberOfEdges(), 1); // 0 - 1
    }

    {
        std::unordered_set<node> nodes = {0, 1};
        auto res = GraphTools::subgraphAndNeighborsFromNodes(G, nodes, true);
        EXPECT_EQ(res.numberOfNodes(), 4);
        EXPECT_EQ(res.numberOfEdges(), 4); // 0-1, 0-2, 1-2, 1-3
    }

    G.addEdge(0, 0);

    for (bool compact : compactOptions) {
        std::unordered_set<node> nodes = {0, 1};
        auto res = GraphTools::subgraphFromNodes(G, nodes.begin(), nodes.end(), compact);
        EXPECT_EQ(res.numberOfNodes(), 2);
        EXPECT_EQ(res.numberOfEdges(), 2); // 0 - 1, 0 - 0
    }
}

TEST_P(GraphToolsGTest, testSubgraphFromNodesDirected) {
    auto G = Graph(4, weighted(), true);

    /**
     *      1
     *   /  |  \
     * 0    |    3
     *   \  |  /
     *      2
     */

    G.addEdge(0, 1, 1.0);
    G.addEdge(0, 2, 2.0);
    G.addEdge(3, 1, 4.0);
    G.addEdge(3, 2, 5.0);
    G.addEdge(1, 2, 3.0);

    std::array<bool, 2> compactOptions{{true, false}};

    for (bool compact : compactOptions) {
        std::unordered_set<node> nodes = {0};
        auto res = GraphTools::subgraphFromNodes(G, nodes.begin(), nodes.end(), compact);

        EXPECT_EQ(weighted(), res.isWeighted());
        EXPECT_TRUE(res.isDirected());

        EXPECT_EQ(res.numberOfNodes(), 1);
        EXPECT_EQ(res.numberOfEdges(), 0);

        if (compact) {
            EXPECT_EQ(res.upperNodeIdBound(), 1);
        } else {
            EXPECT_EQ(res.upperNodeIdBound(), G.upperNodeIdBound());
        }
    }

    {
        std::unordered_set<node> nodes = {0};
        auto res = GraphTools::subgraphAndNeighborsFromNodes(G, nodes, true);
        EXPECT_EQ(res.numberOfNodes(), 3);
        EXPECT_EQ(res.numberOfEdges(), 2); // 0->1, 0->2, NOT 1->2
    }

    for (bool compact : compactOptions) {
        std::unordered_set<node> nodes = {0, 1};
        auto res = GraphTools::subgraphFromNodes(G, nodes.begin(), nodes.end(), compact);
        EXPECT_EQ(res.numberOfNodes(), 2);
        EXPECT_EQ(res.numberOfEdges(), 1); // 0 -> 1
    }

    {
        std::unordered_set<node> nodes = {0, 1};
        auto res = GraphTools::subgraphAndNeighborsFromNodes(G, nodes, true);
        EXPECT_EQ(res.numberOfNodes(), 3);
        EXPECT_EQ(res.numberOfEdges(), 3); // 0->1, 0->2, 1->2
    }

    {
        std::unordered_set<node> nodes = {0, 1};
        auto res = GraphTools::subgraphAndNeighborsFromNodes(G, nodes, true, true);
        EXPECT_EQ(res.numberOfNodes(), 4);
        EXPECT_EQ(res.numberOfEdges(), 4); // 0->1, 0->2, 1->2, 3->1
    }
}

TEST_P(GraphToolsGTest, testTranspose) {
    auto G = Graph(4, weighted(), true);

    /**
     *      1
     *   /  |  \
     * 0    |    3
     *   \  |  /
     *      2
     */

    G.addNode(); // node 4
    G.addNode(); // node 5
    G.addNode(); // node 6
    G.removeNode(5);

    G.addEdge(0, 0, 3.14);
    G.addEdge(0, 4, 3.14);
    G.removeEdge(0, 4);
    G.addEdge(0, 6, 3.14);

    // expect throw error when G is undirected
    if (!G.isDirected()) {
        EXPECT_ANY_THROW(GraphTools::transpose(G));
    } else {
        Graph Gtrans = GraphTools::transpose(G);
        // check summation statistics
        EXPECT_EQ(G.numberOfNodes(), Gtrans.numberOfNodes());
        EXPECT_EQ(G.upperNodeIdBound(), Gtrans.upperNodeIdBound());
        EXPECT_EQ(G.numberOfEdges(), Gtrans.numberOfEdges());
        EXPECT_EQ(G.upperEdgeIdBound(), Gtrans.upperEdgeIdBound());
        EXPECT_EQ(G.totalEdgeWeight(), Gtrans.totalEdgeWeight());
        EXPECT_EQ(G.numberOfSelfLoops(), Gtrans.numberOfSelfLoops());

        // test for regular edges
        EXPECT_TRUE(G.hasEdge(0, 6));
        EXPECT_FALSE(G.hasEdge(6, 0));
        EXPECT_TRUE(Gtrans.hasEdge(6, 0));
        EXPECT_FALSE(Gtrans.hasEdge(0, 6));
        // .. and for selfloops
        EXPECT_TRUE(G.hasEdge(0, 0));
        EXPECT_TRUE(Gtrans.hasEdge(0, 0));

        // check for edge weights
        EXPECT_EQ(G.weight(0, 6), weighted() ? 3.14 : defaultEdgeWeight);
        EXPECT_EQ(Gtrans.weight(6, 0), weighted() ? 3.14 : defaultEdgeWeight);
        EXPECT_EQ(G.weight(0, 0), Gtrans.weight(0, 0));
    }
}

TEST_P(GraphToolsGTest, testToUndirected) {
    constexpr count n = 200;
    constexpr double p = 0.2;

    auto testGraphs = [&](const Graph &G, const Graph &G1) {
        EXPECT_EQ(G.numberOfNodes(), G1.numberOfNodes());
        EXPECT_EQ(G.upperNodeIdBound(), G1.upperNodeIdBound());
        EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        EXPECT_EQ(G.upperEdgeIdBound(), G1.upperEdgeIdBound());
        EXPECT_EQ(G.isWeighted(), G1.isWeighted());
        EXPECT_NE(G.isDirected(), G1.isDirected());
        EXPECT_EQ(G.hasEdgeIds(), G1.hasEdgeIds());

        G.forEdges([&](node u, node v, edgeweight w) {
            EXPECT_TRUE(G1.hasEdge(u, v));
            EXPECT_DOUBLE_EQ(G1.weight(u, v), w);
        });
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(n, p, true).generate();
        if (weighted()) {
            G = generateRandomWeights(G);
        }
        auto G1 = GraphTools::toUndirected(G);
        testGraphs(G, G1);
    }
}

TEST_P(GraphToolsGTest, testToUnWeighted) {
    constexpr count n = 200;
    constexpr double p = 0.2;

    auto testGraphs = [&](const Graph &G, const Graph &G1) {
        EXPECT_EQ(G.numberOfNodes(), G1.numberOfNodes());
        EXPECT_EQ(G.upperNodeIdBound(), G1.upperNodeIdBound());
        EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        EXPECT_NE(G.isWeighted(), G1.isWeighted());
        EXPECT_EQ(G.isDirected(), G1.isDirected());
        EXPECT_EQ(G.hasEdgeIds(), G1.hasEdgeIds());

        G.forEdges([&](node u, node v) {
            EXPECT_TRUE(G1.hasEdge(u, v));
            if (G1.isWeighted()) {
                EXPECT_EQ(G1.weight(u, v), defaultEdgeWeight);
            }
        });
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(n, p, directed()).generate();
        auto G1 = GraphTools::toWeighted(G);
        testGraphs(G, G1);

        G = generateRandomWeights(G);
        G1 = GraphTools::toUnweighted(G);
        testGraphs(G, G1);
    }
}

TEST_P(GraphToolsGTest, testAppend) {
    constexpr count n1 = 100, n2 = 50;
    constexpr double p1 = 0.01, p2 = 0.05;
    constexpr count nodesToDelete = 20;

    auto testGraphs = [&](const Graph &G, const Graph &G1, const Graph &G2) {
        EXPECT_EQ(G.numberOfNodes(), G1.numberOfNodes() + G2.numberOfNodes());
        EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges() + G2.numberOfEdges());
        EXPECT_EQ(G.isDirected(), G1.isDirected());
        EXPECT_EQ(G.isDirected(), G2.isDirected());
        EXPECT_EQ(G.isWeighted(), G1.isWeighted());
        EXPECT_EQ(G.isWeighted(), G2.isWeighted());

        std::unordered_map<node, node> nodeMap;
        node v = G1.upperNodeIdBound();
        for (node u = 0; u < G2.upperNodeIdBound(); ++u) {
            if (G2.hasNode(u)) {
                nodeMap[u] = v++;
            }
        }

        G1.forNodes([&](node u) { EXPECT_TRUE(G.hasNode(u)); });
        G1.forEdges([&](node u, node v) { EXPECT_TRUE(G.hasEdge(u, v)); });
        G2.forNodes([&](node u) { EXPECT_TRUE(G.hasNode(nodeMap[u])); });
        G2.forEdges([&](node u, node v) { EXPECT_TRUE(G.hasEdge(nodeMap[u], nodeMap[v])); });
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G1 = ErdosRenyiGenerator(n1, p1, directed()).generate();
        auto G2 = ErdosRenyiGenerator(n2, p2, directed()).generate();

        if (weighted()) {
            G1 = generateRandomWeights(G1);
            G2 = generateRandomWeights(G2);
        }

        auto G(G1);
        GraphTools::append(G, G2);
        testGraphs(G, G1, G2);

        for (count i = 0; i < nodesToDelete; ++i) {
            G1.removeNode(GraphTools::randomNode(G1));
            G2.removeNode(GraphTools::randomNode(G2));
            auto G3(G1);
            GraphTools::append(G3, G2);
            testGraphs(G3, G1, G2);
        }
    }
}

TEST_P(GraphToolsGTest, testMerge) {
    constexpr count n1 = 100, n2 = 150;
    constexpr double p1 = 0.01, p2 = 0.05;

    auto testGraphs = [&](const Graph &Gorig, const Graph &Gmerge, const Graph &G1) {
        for (node u = 0; u < std::max(Gorig.upperNodeIdBound(), G1.upperNodeIdBound()); ++u) {
            EXPECT_EQ(Gmerge.hasNode(u), Gorig.hasNode(u) || G1.hasNode(u));
        }
        Gorig.forEdges([&](node u, node v) { EXPECT_TRUE(Gmerge.hasEdge(u, v)); });
        G1.forEdges([&](node u, node v) { EXPECT_TRUE(Gmerge.hasEdge(u, v)); });

        Gmerge.forEdges([&](node u, node v, edgeweight w) {
            if (Gorig.hasNode(u) && Gorig.hasNode(v) && Gorig.hasEdge(u, v)) {
                EXPECT_EQ(Gorig.weight(u, v), w);
            } else {
                EXPECT_EQ(G1.weight(u, v), w);
            }
        });
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto Gorig = ErdosRenyiGenerator(n1, p1, directed()).generate();
        auto G1 = ErdosRenyiGenerator(n2, p2, directed()).generate();

        if (weighted()) {
            Gorig = generateRandomWeights(Gorig);
            G1 = generateRandomWeights(G1);
        }

        auto Gmerge(Gorig);
        GraphTools::merge(Gmerge, G1);
        testGraphs(Gorig, Gmerge, G1);
    }
}

TEST_P(GraphToolsGTest, testEdgesSortedByWeight) {
    const auto hasEdgesSortedByWeight = [](const Graph &G, bool decreasing) -> bool {
        for (node u : G.nodeRange()) {
            node prevNode = decreasing ? none : 0;
            edgeweight prevWeight =
                (decreasing ? 1. : -1.) * std::numeric_limits<edgeweight>::max();
            bool sorted = true;
            G.forNeighborsOf(u, [&](node v, edgeweight ew) {
                if (ew == prevWeight)
                    sorted = prevNode <= v;
                else
                    sorted = decreasing ? ew < prevWeight : ew > prevWeight;
                prevNode = v;
                prevWeight = ew;
            });

            if (!sorted)
                return false;
        }

        return true;
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);
        auto G = METISGraphReader{}.read("input/PGPgiantcompo.graph");

        if (weighted()) {
            G = GraphTools::toWeighted(G);
            G.forEdges([&G](node u, node v) { G.setWeight(u, v, Aux::Random::real(10)); });
        }

        GraphTools::sortEdgesByWeight(G);
        EXPECT_TRUE(hasEdgesSortedByWeight(G, false));
        EXPECT_TRUE(G.checkConsistency());
        GraphTools::sortEdgesByWeight(G, true);
        EXPECT_TRUE(hasEdgesSortedByWeight(G, true));
        EXPECT_TRUE(G.checkConsistency());
    }
}
} // namespace NetworKit
