// no-networkit-format
/*
 * CentralityGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <iomanip>
#include <iostream>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/centrality/ApproxBetweenness.hpp>
#include <networkit/centrality/ApproxCloseness.hpp>
#include <networkit/centrality/ApproxElectricalCloseness.hpp>
#include <networkit/centrality/ApproxSpanningEdge.hpp>
#include <networkit/centrality/ApproxGroupBetweenness.hpp>
#include <networkit/centrality/Betweenness.hpp>
#include <networkit/centrality/Closeness.hpp>
#include <networkit/centrality/CoreDecomposition.hpp>
#include <networkit/centrality/DegreeCentrality.hpp>
#include <networkit/centrality/DynApproxBetweenness.hpp>
#include <networkit/centrality/DynKatzCentrality.hpp>
#include <networkit/centrality/DynTopHarmonicCloseness.hpp>
#include <networkit/centrality/EigenvectorCentrality.hpp>
#include <networkit/centrality/EstimateBetweenness.hpp>
#include <networkit/centrality/GedWalk.hpp>
#include <networkit/centrality/GroupCloseness.hpp>
#include <networkit/centrality/GroupClosenessGrowShrink.hpp>
#include <networkit/centrality/GroupClosenessLocalSwaps.hpp>
#include <networkit/centrality/GroupClosenessLocalSearch.hpp>
#include <networkit/centrality/GroupDegree.hpp>
#include <networkit/centrality/GroupHarmonicCloseness.hpp>
#include <networkit/centrality/HarmonicCloseness.hpp>
#include <networkit/centrality/KPathCentrality.hpp>
#include <networkit/centrality/KadabraBetweenness.hpp>
#include <networkit/centrality/KatzCentrality.hpp>
#include <networkit/centrality/LaplacianCentrality.hpp>
#include <networkit/centrality/LocalClusteringCoefficient.hpp>
#include <networkit/centrality/PageRank.hpp>
#include <networkit/centrality/PermanenceCentrality.hpp>
#include <networkit/centrality/SpanningEdgeCentrality.hpp>
#include <networkit/centrality/TopCloseness.hpp>
#include <networkit/centrality/TopHarmonicCloseness.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/Dijkstra.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>
#include <networkit/structures/Cover.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class CentralityGTest : public testing::TestWithParam<std::pair<bool, bool>> {
protected:
    bool isDirected() const noexcept;
    bool isWeighted() const noexcept;
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, CentralityGTest,
        testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                        std::make_pair(false, true),
                        std::make_pair(true, true)));

bool CentralityGTest::isWeighted() const noexcept { return GetParam().first; }

bool CentralityGTest::isDirected() const noexcept { return GetParam().second; }

TEST_F(CentralityGTest, testBetweennessCentrality) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    Betweenness centrality(G);
    centrality.run();
    std::vector<double> bc = centrality.scores();
    std::vector<double> result = {0., 0., 15., 3., 3., 1.};

    const double tol = 1e-3;
    G.forNodes([&](node u) {
        EXPECT_NEAR(result[u], bc[u], tol);
    });

    Betweenness centralityNorm(G, true);
    centralityNorm.run();
    bc = centralityNorm.scores();
    const double pairs = (n - 1) * (n - 2);
    G.forNodes([&](node u) {
        EXPECT_NEAR(result[u] / pairs, bc[u], tol);
    });
}

TEST_F(CentralityGTest, testBetweenness2Centrality) {
    /* Graph:
          0    3
          \  / \
          2    5
          /  \ /
          1    4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    Betweenness centrality(G);
    centrality.run();
    std::vector<double> bc = centrality.scores();

    const double tol = 1e-3;
    EXPECT_NEAR(0.0, bc[0], tol);
    EXPECT_NEAR(0.0, bc[1], tol);
    EXPECT_NEAR(15.0, bc[2], tol);
    EXPECT_NEAR(3.0, bc[3], tol);
    EXPECT_NEAR(3.0, bc[4], tol);
    EXPECT_NEAR(1.0, bc[5], tol);
}

TEST_F(CentralityGTest, runApproxBetweennessSmallGraph) {
    /* Graph:
     0    3
      \  / \
       2   5
      / \ /
     1   4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    double epsilon = 0.1; // error
    double delta = 0.1;   // confidence
    ApproxBetweenness centrality(G, epsilon, delta);
    centrality.run();
    std::vector<double> bc = centrality.scores();

    ASSERT_LE(centrality.scores().size(), 1.0 / (epsilon * epsilon));

    DEBUG("scores: ", bc);
}

TEST_F(CentralityGTest, testBetweennessCentralityWeighted) {
    /* Graph:
     0    3   6
      \  / \ /
       2 -- 5
      /  \ / \
     1    4   7

     Edges in the upper row have weight 3,
     the edge in the middle row has weight 1.5,
     edges in the lower row have weight 2.
    */
    count n = 8;
    Graph G(n, true);

    G.addEdge(0, 2, 3);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(2, 4, 2);
    G.addEdge(2, 5, 1.5);
    G.addEdge(3, 5, 3);
    G.addEdge(4, 5, 2);
    G.addEdge(5, 6, 3);
    G.addEdge(5, 7, 2);

    Betweenness centrality(G);
    centrality.run();
    std::vector<double> bc = centrality.scores();

    const double tol = 1e-3;
    EXPECT_NEAR(0.0, bc[0], tol);
    EXPECT_NEAR(0.0, bc[1], tol);
    EXPECT_NEAR(23.0, bc[2], tol);
    EXPECT_NEAR(0.0, bc[3], tol);
    EXPECT_NEAR(0.0, bc[4], tol);
    EXPECT_NEAR(23.0, bc[5], tol);
    EXPECT_NEAR(0.0, bc[6], tol);
    EXPECT_NEAR(0.0, bc[7], tol);
}

// TODO: replace by smaller graph
TEST_F(CentralityGTest, testKatzCentralityDirected) {
    const auto G = SNAPGraphReader{}.read("input/wiki-Vote.txt");
    KatzCentrality kc(G, 5e-4);

    DEBUG("start kc run");
    kc.run();
    DEBUG("finish kc");

    EXPECT_EQ(kc.ranking().front().first, 699);
}

TEST_F(CentralityGTest, testKatzTopk) {
    const auto G = METISGraphReader{}.read("input/caidaRouterLevel.graph");

    KatzCentrality exactAlgo(G, 0, 1.0);
    DynKatzCentrality topAlgo(G, 100);
    exactAlgo.run();
    topAlgo.run();

    // We cannot compare the ranking as the algorithms might return different
    // rankings for nodes that have equal/nearly equal scores. Instead,
    // epsilon-compare the exact scores of the i-th node and the expected i-th
    // node.
    auto exactRanking = exactAlgo.ranking();
    auto topRanking = topAlgo.ranking();
    for (count i = 0; i < std::min(G.numberOfNodes(), static_cast<count>(100)); i++)
        EXPECT_NEAR(exactAlgo.score(topRanking[i].first),
                    exactRanking[i].second, 1e-6);
}

TEST_F(CentralityGTest, testKatzDynamicAddition) {
    Graph G = METISGraphReader{}.read("input/caidaRouterLevel.graph");
    DynKatzCentrality kc(G, 100);
    DEBUG("start kc run");
    kc.run();
    DEBUG("finish kc");
    Aux::Random::setSeed(42, false);
    node u, v;
    do {
        u = GraphTools::randomNode(G);
        v = GraphTools::randomNode(G);
    } while (G.hasEdge(u, v));
    GraphEvent e(GraphEvent::EDGE_ADDITION, u, v, 1.0);
    kc.update(e);
    G.addEdge(u, v);
    DynKatzCentrality kc2(G, 100);
    kc2.run();
    const edgeweight tol = 1e-9;
    for (count i = 0; i <= std::min(kc.levelReached, kc2.levelReached); i++) {
        INFO("i = ", i);
        G.forNodes([&](node u) {
            EXPECT_EQ(kc.nPaths[i][u], kc2.nPaths[i][u]);
        });
    }
    G.forNodes([&](node u) {
        EXPECT_NEAR(kc.score(u), kc2.score(u), tol);
        EXPECT_NEAR(kc.bound(u), kc2.bound(u), tol);
    });

    INFO("Level reached: ", kc.levelReached, ", ", kc2.levelReached);
}

TEST_F(CentralityGTest, testKatzDynamicDeletion) {
    Graph G = METISGraphReader{}.read("input/caidaRouterLevel.graph");
    DynKatzCentrality kc(G, 100);
    DEBUG("start kc run");
    kc.run();
    DEBUG("finish kc");
    std::pair<node, node> p = GraphTools::randomEdge(G);
    node u = p.first;
    node v = p.second;
    INFO("Deleting edge ", u, ", ", v);
    GraphEvent e(GraphEvent::EDGE_REMOVAL, u, v, 1.0);
    G.removeEdge(u, v);
    kc.update(e);
    DynKatzCentrality kc2(G, 100);
    kc2.run();
    const edgeweight tol = 1e-9;
    for (count i = 0; i <= std::min(kc.levelReached, kc2.levelReached); i++) {
        INFO("i = ", i);
        G.forNodes([&](node u) {
            if (kc.nPaths[i][u] != kc2.nPaths[i][u]) {
                INFO("i = ", i, ", node ", u,
                     ", dyn kc paths: ", kc.nPaths[i][u],
                     ", stat paths: ", kc2.nPaths[i][u]);
            }
            EXPECT_EQ(kc.nPaths[i][u], kc2.nPaths[i][u]);
        });
    }
    G.forNodes([&](node u) {
        EXPECT_NEAR(kc.score(u), kc2.score(u), tol);
        EXPECT_NEAR(kc.bound(u), kc2.bound(u), tol);
    });

    INFO("Level reached: ", kc.levelReached, ", ", kc2.levelReached);
}

TEST_F(CentralityGTest, testKatzDynamicBuilding) {
    Graph GIn = METISGraphReader{}.read("input/hep-th.graph");

    // Find a single max-degree node and add its edges to G.
    // (This guarantees that alpha is correct.)
    node maxNode = 0;
    GIn.forNodes([&](node u) {
        if (GIn.degree(u) > GIn.degree(maxNode))
            maxNode = u;
    });

    Graph G(GIn.upperNodeIdBound());

    GIn.forEdgesOf(maxNode, [&](node u, edgeweight) { G.addEdge(maxNode, u); });

    // Now run the algo. and add other some edges to check the correctness of
    // the dynamic part.
    DynKatzCentrality dynAlgo(G, 100);
    dynAlgo.run();

    count edgesProcessed = 0;
    GIn.forEdges([&](node u, node v) {
        if (u == maxNode || v == maxNode)
            return;
        if (edgesProcessed > 1000)
            return;
        GraphEvent e(GraphEvent::EDGE_ADDITION, u, v, 1.0);
        dynAlgo.update(e);
        G.addEdge(u, v);
        edgesProcessed++;
    });

    DynKatzCentrality topAlgo(G, 100);
    topAlgo.run();

    auto topRanking = topAlgo.ranking();
    auto dynRanking = dynAlgo.ranking();
    for (count i = 0; i < std::min(G.numberOfNodes(), count{100}); i++)
        EXPECT_FALSE(
            dynAlgo.areDistinguished(topRanking[i].first, dynRanking[i].first))
            << "Nodes " << topRanking[i].first << " and " << dynRanking[i].first
            << " should not be distinguished!";
}

TEST_F(CentralityGTest, testKatzDirectedAddition) {
    // Same graph as in testCoreDecompositionDirected.
    count n = 16;
    Graph G(n, false, true);

    G.addEdge(2, 4);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 7);
    G.addEdge(6, 7);

    G.addEdge(6, 8);
    G.addEdge(6, 9);
    G.addEdge(6, 11);
    G.addEdge(7, 12);
    G.addEdge(8, 9);

    G.addEdge(8, 10);
    G.addEdge(8, 11);
    G.addEdge(8, 13);
    G.addEdge(9, 10);
    G.addEdge(9, 11);

    G.addEdge(9, 13);
    G.addEdge(10, 11);
    G.addEdge(10, 13);
    G.addEdge(10, 14);
    G.addEdge(11, 13);

    G.addEdge(11, 14);
    G.addEdge(12, 15);
    G.addEdge(13, 14);
    G.addEdge(14, 15);

    DynKatzCentrality kc(G, 5);
    kc.run();

    node u, v;
    Aux::Random::setSeed(42, false);
    do {
        u = GraphTools::randomNode(G);
        v = GraphTools::randomNode(G);
    } while (G.hasEdge(u, v));
    GraphEvent e(GraphEvent::EDGE_ADDITION, u, v, 1.0);
    kc.update(e);
    G.addEdge(u, v);

    DynKatzCentrality kc2(G, 5);
    kc2.run();

    for (count i = 0; i <= std::min(kc.levelReached, kc2.levelReached); i++) {
        G.forNodes(
            [&](node u) { EXPECT_EQ(kc.nPaths[i][u], kc2.nPaths[i][u]); });
    }
    const edgeweight tol = 1e-9;
    G.forNodes([&](node u) {
        EXPECT_NEAR(kc.score(u), kc2.score(u), tol);
        EXPECT_NEAR(kc.bound(u), kc2.bound(u), tol);
    });

    INFO("Level reached: ", kc.levelReached, ", ", kc2.levelReached);
}

TEST_F(CentralityGTest, testKatzDirectedDeletion) {
    // Same graph as in testCoreDecompositionDirected.
    count n = 16;
    Graph G(n, false, true);

    G.addEdge(2, 4);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 7);
    G.addEdge(6, 7);

    G.addEdge(6, 8);
    G.addEdge(6, 9);
    G.addEdge(6, 11);
    G.addEdge(7, 12);
    G.addEdge(8, 9);

    G.addEdge(8, 10);
    G.addEdge(8, 11);
    G.addEdge(8, 13);
    G.addEdge(9, 10);
    G.addEdge(9, 11);

    G.addEdge(9, 13);
    G.addEdge(10, 11);
    G.addEdge(10, 13);
    G.addEdge(10, 14);
    G.addEdge(11, 13);

    G.addEdge(11, 14);
    G.addEdge(12, 15);
    G.addEdge(13, 14);
    G.addEdge(14, 15);

    DynKatzCentrality kc(G, 5);
    kc.run();

    Aux::Random::setSeed(42, false);
    std::pair<node, node> p = GraphTools::randomEdge(G);
    node u = p.first;
    node v = p.second;
    INFO("Removing ", u, " -> ", v);
    GraphEvent e(GraphEvent::EDGE_REMOVAL, u, v, 1.0);
    G.removeEdge(u, v);
    kc.update(e);

    DynKatzCentrality kc2(G, 5);
    kc2.run();

    for (count i = 0; i <= std::min(kc.levelReached, kc2.levelReached); i++) {
        G.forNodes([&](node u) {
            EXPECT_EQ(kc.nPaths[i][u], kc2.nPaths[i][u])
                << i << "-length paths ending in " << u << " do not match!";
        });
    }
    const edgeweight tol = 1e-9;
    G.forNodes([&](node u) {
        EXPECT_NEAR(kc.score(u), kc2.score(u), tol);
        EXPECT_NEAR(kc.bound(u), kc2.bound(u), tol);
    });

    INFO("Level reached: ", kc.levelReached, ", ", kc2.levelReached);
}

TEST_P(CentralityGTest, testPageRank) {
    SNAPGraphReader reader(isDirected());
    auto G = reader.read("input/wiki-Vote.txt");

    auto doTest = [&G](PageRank::Norm norm) {
        PageRank pr(G);
        pr.norm = norm;
        pr.run();

        auto pr_ranking = pr.ranking();
        constexpr double eps = 1e-3;
        if (G.isDirected()) {
            EXPECT_EQ(pr_ranking[0].first, 326);
            EXPECT_NEAR(pr_ranking[0].second, 0.00460, eps);
        } else {
            EXPECT_EQ(pr_ranking[0].first, 699);
            EXPECT_NEAR(pr_ranking[0].second, 0.00432, eps);
        }

        constexpr count maxIterations = 2;
        pr.maxIterations = maxIterations;
        pr.run();
        EXPECT_LE(pr.numberOfIterations(), maxIterations);
    };

    doTest(PageRank::Norm::L1Norm);
    doTest(PageRank::Norm::L2Norm);
}

TEST_P(CentralityGTest, testNormalizedPageRank) {
    SNAPGraphReader reader(isDirected());
    auto G = reader.read("input/wiki-Vote.txt");

    auto doTest = [&G](PageRank::Norm norm) {
        PageRank pr(G, 0.85, 1e-8, true);
        pr.norm = norm;
        pr.run();

        auto pr_ranking = pr.ranking();
        const double tol = 1e-4;
        if (G.isDirected()) {
            EXPECT_EQ(pr_ranking[0].first, 326);
            EXPECT_NEAR(pr_ranking[0].second, 57.6657, tol);
        } else {
            EXPECT_EQ(pr_ranking[0].first, 699);
            EXPECT_NEAR(pr_ranking[0].second, 205.7324, tol);
        }
        const count maxIterations = 2;
        pr.maxIterations = maxIterations;
        pr.run();
        EXPECT_LE(pr.numberOfIterations(), maxIterations);
    };

    doTest(PageRank::Norm::L1Norm);
    doTest(PageRank::Norm::L2Norm);
}

TEST_F(CentralityGTest, testEigenvectorCentrality) {
    /* Graph:
     0    3   6
      \  / \ /
       2 -- 5
      /  \ / \
     1    4   7

     Edges in the upper row have weight 3,
     the edge in the middle row has weight 1.5,
     edges in the lower row have weight 2.
    */
    count n = 8;
    Graph G(n, true);

    G.addEdge(0, 2, 3);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(2, 4, 2);
    G.addEdge(2, 5, 1.5);
    G.addEdge(3, 5, 3);
    G.addEdge(4, 5, 2);
    G.addEdge(5, 6, 3);
    G.addEdge(5, 7, 2);

    EigenvectorCentrality centrality(G);
    centrality.run();
    std::vector<double> cen = centrality.scores();

    // computed with Matlab
    const double tol = 1e-4;
    EXPECT_NEAR(0.2254, fabs(cen[0]), tol);
    EXPECT_NEAR(0.1503, fabs(cen[1]), tol);
    EXPECT_NEAR(0.5290, fabs(cen[2]), tol);
    EXPECT_NEAR(0.4508, fabs(cen[3]), tol);
    EXPECT_NEAR(0.3006, fabs(cen[4]), tol);
    EXPECT_NEAR(0.5290, fabs(cen[5]), tol);
    EXPECT_NEAR(0.2254, fabs(cen[6]), tol);
    EXPECT_NEAR(0.1503, fabs(cen[7]), tol);
}

TEST_F(CentralityGTest, testPageRankCentrality) {
    /* Graph:
     0    3   6
      \  / \ /
       2 -- 5
      /  \ / \
     1    4   7

     Edges in the upper row have weight 3,
     the edge in the middle row has weight 1.5,
     edges in the lower row have weight 2.
    */
    count n = 8;
    Graph G(n, true);

    G.addEdge(0, 2, 3);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(2, 4, 2);
    G.addEdge(2, 5, 1.5);
    G.addEdge(3, 5, 3);
    G.addEdge(4, 5, 2);
    G.addEdge(5, 6, 3);
    G.addEdge(5, 7, 2);

    double damp = 0.85;
    PageRank centrality(G, damp);
    centrality.run();
    std::vector<double> cen = centrality.scores();

    // compare to Matlab results
    const double tol = 1e-4;
    EXPECT_NEAR(0.0753, fabs(cen[0]), tol);
    EXPECT_NEAR(0.0565, fabs(cen[1]), tol);
    EXPECT_NEAR(0.2552, fabs(cen[2]), tol);
    EXPECT_NEAR(0.1319, fabs(cen[3]), tol);
    EXPECT_NEAR(0.0942, fabs(cen[4]), tol);
    EXPECT_NEAR(0.2552, fabs(cen[5]), tol);
    EXPECT_NEAR(0.0753, fabs(cen[6]), tol);
    EXPECT_NEAR(0.0565, fabs(cen[7]), tol);
}

TEST_F(CentralityGTest, testNormalizedPageRankCentrality) {
    /* Graph:
     0    3   6
      \  / \ /
       2 -- 5
      /  \ / \
     1    4   7

     Edges in the upper row have weight 3,
     the edge in the middle row has weight 1.5,
     edges in the lower row have weight 2.
    */
    count n = 8;
    Graph G(n, true);

    G.addEdge(0, 2, 3);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(2, 4, 2);
    G.addEdge(2, 5, 1.5);
    G.addEdge(3, 5, 3);
    G.addEdge(4, 5, 2);
    G.addEdge(5, 6, 3);
    G.addEdge(5, 7, 2);

    double damp = 0.85;
    PageRank centrality(G, damp, 1e-8, true);
    centrality.run();
    std::vector<double> cen = centrality.scores();

    // compare to Matlab results
    const double tol = 1e-4;
    EXPECT_NEAR(4.0175, fabs(cen[0]), tol);
    EXPECT_NEAR(3.0117, fabs(cen[1]), tol);
    EXPECT_NEAR(13.6083, fabs(cen[2]), tol);
    EXPECT_NEAR(7.0350, fabs(cen[3]), tol);
    EXPECT_NEAR(5.0233, fabs(cen[4]), tol);
    EXPECT_NEAR(13.6083, fabs(cen[5]), tol);
    EXPECT_NEAR(4.0175, fabs(cen[6]), tol);
    EXPECT_NEAR(3.0117, fabs(cen[7]), tol);
}

TEST_F(CentralityGTest, benchSequentialBetweennessCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    Betweenness bc(G);
    bc.run();
    std::vector<std::pair<node, double>> ranking = bc.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(CentralityGTest, benchParallelBetweennessCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    Betweenness bc(G);
    bc.run();
    std::vector<std::pair<node, double>> ranking = bc.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(CentralityGTest, benchEigenvectorCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    EigenvectorCentrality cen(G);
    cen.run();
    std::vector<std::pair<node, double>> ranking = cen.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(CentralityGTest, benchPageRankCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    double damp = 0.85;
    PageRank cen(G, damp);
    cen.run();
    std::vector<std::pair<node, double>> ranking = cen.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(CentralityGTest, benchNormalizedPageRankCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    double damp = 0.85;
    PageRank cen(G, damp, 1e-8, true);
    cen.run();
    std::vector<std::pair<node, double>> ranking = cen.ranking();
    INFO("Highest rank: ", ranking[0].first, " with score ", ranking[0].second);
}

TEST_F(CentralityGTest, runEstimateBetweenness) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");

    EstimateBetweenness abc2(G, 100);
    abc2.run();

    DEBUG("approximated betweenness scores: ", abc2.scores());
}

// FIXME look out for tolerance limit in paper sample nodes
TEST_F(CentralityGTest, testApproxClosenessCentralityOnRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");

    ApproxCloseness acc(G, 453, 0, true);
    acc.run();

    std::vector<double> acc_scores = acc.scores();

    ASSERT_EQ(acc_scores.size(), 453);

    // compare sampled approx closeness values vs real closeness values
    // its a good compromise to not let the exact closness algorithm run
    ASSERT_NEAR(acc_scores[0], 0.416206, 0.000001);
    ASSERT_NEAR(acc_scores[10], 0.355906, 0.000001);
    ASSERT_NEAR(acc_scores[87], 0.420465, 0.000001);
    ASSERT_NEAR(acc_scores[121], 0.38865, 0.000001);
    ASSERT_NEAR(acc_scores[178], 0.4, 0.000001);
    ASSERT_NEAR(acc_scores[254], 0.397188, 0.000001);
    ASSERT_NEAR(acc_scores[307], 0.398238, 0.000001);
    ASSERT_NEAR(acc_scores[398], 0.37604, 0.000001);
    ASSERT_NEAR(acc_scores[406], 0.360734, 0.000001);
    ASSERT_NEAR(acc_scores[446], 0.396491, 0.000001);
}

TEST_F(CentralityGTest, testApproxClosenessCentralityOnToyGraph) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    ApproxCloseness acc(G, 6, 0.1, false);
    acc.run();
    std::vector<double> cc = acc.scores();

    double maximum = acc.maximum();

    const double tol = 0.2;
    EXPECT_NEAR(0.1, cc[0], tol);
    EXPECT_NEAR(0.1, cc[1], tol);
    EXPECT_NEAR(0.166667, cc[2], tol);
    EXPECT_NEAR(0.125, cc[3], tol);
    EXPECT_NEAR(0.125, cc[4], tol);
    EXPECT_NEAR(0.1, cc[5], tol);
    EXPECT_NEAR(0.2, maximum, tol);

    ApproxCloseness acc2(G, 4, 0.1, true);
    acc2.run();
    std::vector<double> cc2 = acc2.scores();

    double maximum2 = acc2.maximum();

    EXPECT_NEAR(0.5, cc2[0], tol);
    EXPECT_NEAR(0.5, cc2[1], tol);
    EXPECT_NEAR(0.833335, cc2[2], tol);
    EXPECT_NEAR(0.625, cc2[3], tol);
    EXPECT_NEAR(0.625, cc2[4], tol);
    EXPECT_NEAR(0.5, cc2[5], tol);
    EXPECT_NEAR(0.2, maximum2, tol);
}

TEST_F(CentralityGTest, testEdgeBetweennessCentrality) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.indexEdges();

    Betweenness centrality(G, false, true);
    centrality.run();
    std::vector<double> bc = centrality.edgeScores();

    const double tol = 1e-3;
    EXPECT_NEAR(10.0, bc[0], tol);
    EXPECT_NEAR(10.0, bc[1], tol);
    EXPECT_NEAR(10.0, bc[2], tol);
    EXPECT_NEAR(10.0, bc[3], tol);
    EXPECT_NEAR(6.0, bc[4], tol);
    EXPECT_NEAR(6.0, bc[5], tol);
}

TEST_F(CentralityGTest, debugEdgeBetweennessCentrality) {
    auto path = "input/PGPgiantcompo.graph";
    METISGraphReader reader;
    Graph G = reader.read(path);
    G.indexEdges();

    Betweenness centrality(G, false, true);
    centrality.run();
    std::vector<double> bc = centrality.edgeScores();
}

TEST_P(CentralityGTest, testClosenessCentrality) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n, isWeighted(), isDirected());

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    if (isDirected()) {
        G.addEdge(4, 1);
        G.addEdge(3, 0);
        G.addEdge(5, 2);
    }

    if (isWeighted()) {
        Aux::Random::setSeed(42, false);
        G.forEdges([&](node u, node v) { G.setWeight(u, v, Aux::Random::probability()); });
    }

    auto computeCloseness = [&](node u, ClosenessVariant variant, bool normalized) {
        Dijkstra dij(G, u, false);
        dij.run();
        auto dists = dij.getDistances();

        double sumDist = 0.;
        double reached = 0.;

        G.forNodes([&](node v) {
            if (variant == ClosenessVariant::standard) {
                sumDist += dists[v];
            } else if (dists[v] != std::numeric_limits<double>::max()) {
                ++reached;
                sumDist += dists[v];
            }
        });

        double score = sumDist;
        if (score) {
            score = (variant == ClosenessVariant::standard) ?
                1.0 / score : (reached - 1.0) / sumDist / (G.numberOfNodes() - 1.0);
        }

        if (normalized) {
            score *= (variant == ClosenessVariant::standard ? G.numberOfNodes() : reached) - 1.0;
        }

        return score;
    };

    for (auto variant : {ClosenessVariant::standard, ClosenessVariant::generalized}) {
        for (auto normalized : {true, false}) {
            Closeness centrality(G, normalized, variant);
            centrality.run();
            const auto bc = centrality.scores();

            G.forNodes([&](node u) {
                EXPECT_DOUBLE_EQ(bc[u], computeCloseness(u, variant, normalized));
            });
        }
    }
}

TEST_F(CentralityGTest, testHarmonicClosenessCentrality) {
    /* Graph:
     0    3
      \  / \
       2    5
      /  \ /
     1    4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    HarmonicCloseness centrality(G, false);
    centrality.run();
    std::vector<double> hc = centrality.scores();

    double maximum = centrality.maximum();

    const double tol = 1e-3;
    EXPECT_NEAR(2.833, hc[0], tol);
    EXPECT_NEAR(2.833, hc[1], tol);
    EXPECT_NEAR(4.5, hc[2], tol);
    EXPECT_NEAR(3.5, hc[3], tol);
    EXPECT_NEAR(3.5, hc[4], tol);
    EXPECT_NEAR(3.1667, hc[5], tol);
    EXPECT_NEAR(1, maximum, tol);
}

TEST_F(CentralityGTest, runKPathCentrality) {
    METISGraphReader reader;
    Graph G = reader.read("input/lesmis.graph");

    KPathCentrality centrality(G);
    centrality.run();
}

TEST_F(CentralityGTest, testCoreDecompositionSimple) {
    count n = 3;
    Graph G(n);
    G.addEdge(0, 1);

    CoreDecomposition coreDec(G);
    coreDec.run();
    std::vector<double> coreness = coreDec.scores();

    EXPECT_EQ(1u, coreness[0]) << "expected coreness";
    EXPECT_EQ(1u, coreness[1]) << "expected coreness";
    EXPECT_EQ(0u, coreness[2]) << "expected coreness";
}

TEST_F(CentralityGTest, testCoreDecomposition) {
    count n = 16;
    Graph G(n);

    // 	// create graph used in Baur et al. and network analysis lecture
    G.addEdge(2, 4);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 7);
    G.addEdge(6, 7);

    G.addEdge(6, 8);
    G.addEdge(6, 9);
    G.addEdge(6, 11);
    G.addEdge(7, 12);
    G.addEdge(8, 9);

    G.addEdge(8, 10);
    G.addEdge(8, 11);
    G.addEdge(8, 13);
    G.addEdge(9, 10);
    G.addEdge(9, 11);

    G.addEdge(9, 13);
    G.addEdge(10, 11);
    G.addEdge(10, 13);
    G.addEdge(10, 14);
    G.addEdge(11, 13);

    G.addEdge(11, 14);
    G.addEdge(12, 15);
    G.addEdge(13, 14);
    G.addEdge(14, 15);

    EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
    EXPECT_EQ(24u, G.numberOfEdges()) << "should have 24 edges";

    // compute core decomposition
    CoreDecomposition coreDec(G);
    coreDec.run();
    std::vector<double> coreness = coreDec.scores();
    // init cores
    // init shells
    Cover cores = coreDec.getCover();
    Partition shells = coreDec.getPartition();

    EXPECT_EQ(0u, coreness[0]) << "expected coreness";
    EXPECT_EQ(0u, coreness[1]) << "expected coreness";
    EXPECT_EQ(1u, coreness[2]) << "expected coreness";
    EXPECT_EQ(1u, coreness[3]) << "expected coreness";
    EXPECT_EQ(1u, coreness[4]) << "expected coreness";
    EXPECT_EQ(1u, coreness[5]) << "expected coreness";
    EXPECT_EQ(3u, coreness[6]) << "expected coreness";
    EXPECT_EQ(2u, coreness[7]) << "expected coreness";
    EXPECT_EQ(4u, coreness[8]) << "expected coreness";
    EXPECT_EQ(4u, coreness[9]) << "expected coreness";
    EXPECT_EQ(4u, coreness[10]) << "expected coreness";
    EXPECT_EQ(4u, coreness[11]) << "expected coreness";
    EXPECT_EQ(2u, coreness[12]) << "expected coreness";
    EXPECT_EQ(4u, coreness[13]) << "expected coreness";
    EXPECT_EQ(3u, coreness[14]) << "expected coreness";
    EXPECT_EQ(2u, coreness[15]) << "expected coreness";

    // for (index e = 0; e < n; e++) {
    // 	EXPECT_EQ(cores.contains(e), true);
    // 	EXPECT_EQ(shells.contains(e), true);
    // }
    // EXPECT_EQ(cores.get, coreness[15]) << "expected coreness";

    // test throw runtime error for self-loop in graph
    Graph H(2);
    H.addEdge(0, 1);
    H.addEdge(1, 1);
    EXPECT_ANY_THROW(CoreDecomposition CoreDec(H));
}

TEST_F(CentralityGTest, benchCoreDecompositionLocal) {
    METISGraphReader reader;
    std::vector<std::string> filenames = {"caidaRouterLevel", "wing",
                                          "astro-ph", "PGPgiantcompo"};

    for (auto f : filenames) {
        std::string filename("input/" + f + ".graph");
        DEBUG("about to read file ", filename);
        Graph G = reader.read(filename);
        G.removeSelfLoops();
        CoreDecomposition coreDec(G, false);
        Aux::Timer timer;
        timer.start();
        coreDec.run();
        timer.stop();
        INFO("Time for ParK of ", filename, ": ", timer.elapsedTag());

        CoreDecomposition coreDec2(G, true);
        timer.start();
        coreDec2.run();
        timer.stop();
        INFO("Time for bucket queue based k-core decomposition of ", filename,
             ": ", timer.elapsedTag());

        G.forNodes(
            [&](node u) { EXPECT_EQ(coreDec.score(u), coreDec2.score(u)); });
    }
}

TEST_F(CentralityGTest, testCoreDecompositionDirected) {
    count n = 16;
    Graph G(n, false, true);

    // 	// create graph used in Baur et al. and network analysis lecture
    G.addEdge(2, 4);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 7);
    G.addEdge(6, 7);

    G.addEdge(6, 8);
    G.addEdge(6, 9);
    G.addEdge(6, 11);
    G.addEdge(7, 12);
    G.addEdge(8, 9);

    G.addEdge(8, 10);
    G.addEdge(8, 11);
    G.addEdge(8, 13);
    G.addEdge(9, 10);
    G.addEdge(9, 11);

    G.addEdge(9, 13);
    G.addEdge(10, 11);
    G.addEdge(10, 13);
    G.addEdge(10, 14);
    G.addEdge(11, 13);

    G.addEdge(11, 14);
    G.addEdge(12, 15);
    G.addEdge(13, 14);
    G.addEdge(14, 15);

    EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
    EXPECT_EQ(24u, G.numberOfEdges()) << "should have 24 edges";

    // compute core decomposition
    CoreDecomposition coreDec(G);
    coreDec.run();
    std::vector<double> coreness = coreDec.scores();

    EXPECT_EQ(0u, coreness[0]) << "expected coreness";
    EXPECT_EQ(0u, coreness[1]) << "expected coreness";
    EXPECT_EQ(1u, coreness[2]) << "expected coreness";
    EXPECT_EQ(1u, coreness[3]) << "expected coreness";
    EXPECT_EQ(1u, coreness[4]) << "expected coreness";
    EXPECT_EQ(1u, coreness[5]) << "expected coreness";
    EXPECT_EQ(3u, coreness[6]) << "expected coreness";
    EXPECT_EQ(2u, coreness[7]) << "expected coreness";
    EXPECT_EQ(4u, coreness[8]) << "expected coreness";
    EXPECT_EQ(4u, coreness[9]) << "expected coreness";
    EXPECT_EQ(4u, coreness[10]) << "expected coreness";
    EXPECT_EQ(4u, coreness[11]) << "expected coreness";
    EXPECT_EQ(2u, coreness[12]) << "expected coreness";
    EXPECT_EQ(4u, coreness[13]) << "expected coreness";
    EXPECT_EQ(3u, coreness[14]) << "expected coreness";
    EXPECT_EQ(2u, coreness[15]) << "expected coreness";
}

TEST_F(CentralityGTest, testLocalClusteringCoefficientUndirected) {
    count n = 16;
    Graph G(n, false, false);

    // 	// create graph used in Baur et al. and network analysis lecture
    G.addEdge(2, 4);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(5, 7);
    G.addEdge(6, 7);

    G.addEdge(6, 8);
    G.addEdge(6, 9);
    G.addEdge(6, 11);
    G.addEdge(7, 12);
    G.addEdge(8, 9);

    G.addEdge(8, 10);
    G.addEdge(8, 11);
    G.addEdge(8, 13);
    G.addEdge(9, 10);
    G.addEdge(9, 11);

    G.addEdge(9, 13);
    G.addEdge(10, 11);
    G.addEdge(10, 13);
    G.addEdge(10, 14);
    G.addEdge(11, 13);

    G.addEdge(11, 14);
    G.addEdge(12, 15);
    G.addEdge(13, 14);
    G.addEdge(14, 15);

    EXPECT_EQ(n, G.numberOfNodes()) << "should have " << n << " vertices";
    EXPECT_EQ(24u, G.numberOfEdges()) << "should have 24 edges";

    // compute core decomposition
    LocalClusteringCoefficient lcc(G);
    lcc.run();
    std::vector<double> lccScores = lcc.scores();
    std::vector<double> reference = {0.0, 0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.5, 0.0,
                                     0.8, 0.8, 0.8, 0.6666666666666666,
                                     0.0, 0.8, 0.5, 0.0};

    EXPECT_EQ(reference, lccScores);

    LocalClusteringCoefficient lccTurbo(G, true);
    lccTurbo.run();
    EXPECT_EQ(reference, lccTurbo.scores());

    // test throw runtime error for self-loop in graph
    Graph H(2);
    H.addEdge(0, 1);
    H.addEdge(1, 1);
    EXPECT_ANY_THROW(LocalClusteringCoefficient lcc(H));
}

TEST_F(CentralityGTest, testLocalClusteringCoefficientUndirected2) {
    Graph G(6, false, false);
    G.addEdge(1, 0);
    G.addEdge(2, 0);
    G.addEdge(2, 1);
    G.addEdge(3, 2);
    G.addEdge(3, 0);
    G.addEdge(3, 1);
    G.addEdge(4, 2);
    G.addEdge(4, 0);
    G.addEdge(5, 3);
    G.addEdge(5, 4);
    G.addEdge(5, 1);
    LocalClusteringCoefficient lcc(G);
    lcc.run();
    std::vector<double> lccScores = lcc.scores();
    std::vector<double> reference = {0.6666666666666666, 0.6666666666666666,
                                     0.6666666666666666, 0.6666666666666666,
                                     0.3333333333333333, 0.3333333333333333};

    EXPECT_EQ(reference, lccScores);
}

TEST_F(CentralityGTest, testSimplePermanence) {
    Graph G(15, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 0);
    G.addEdge(2, 3);
    node v = 4;
    node u = 5;
    G.addEdge(v, 0);
    G.addEdge(v, 1);
    G.addEdge(v, 2);
    G.addEdge(u, 3);
    G.addEdge(u, 2);
    G.addEdge(u, 0);
    G.addEdge(6, 7);
    G.addEdge(7, 8);
    G.addEdge(u, 6);
    G.addEdge(u, 7);
    G.addEdge(u, 8);
    G.addEdge(v, 6);
    G.addEdge(v, 7);
    G.addEdge(9, 10);
    G.addEdge(10, 11);
    G.addEdge(u, 9);
    G.addEdge(v, 10);
    G.addEdge(v, 11);
    G.addEdge(12, 13);
    G.addEdge(13, 14);
    G.addEdge(12, 14);
    G.addEdge(v, 12);
    G.addEdge(v, 14);

    Partition P(G.upperNodeIdBound());
    P.setUpperBound(4);
    P[0] = 0;
    P[1] = 0;
    P[2] = 0;
    P[3] = 0;
    P[v] = 0;
    P[u] = 0;
    P[6] = 1;
    P[7] = 1;
    P[8] = 1;
    P[9] = 2;
    P[10] = 2;
    P[11] = 2;
    P[12] = 3;
    P[13] = 3;
    P[14] = 3;

    ASSERT_EQ(9u, G.degree(v));
    ASSERT_EQ(7u, G.degree(u));

    PermanenceCentrality perm(G, P);
    perm.run();
    EXPECT_DOUBLE_EQ(2.0 / 3.0, perm.getIntraClustering(u));
    EXPECT_DOUBLE_EQ(1, perm.getIntraClustering(v));

    EXPECT_NEAR(-0.19048, perm.getPermanence(u), 0.0005);
    EXPECT_NEAR(0.167, perm.getPermanence(v), 0.0005);
}

TEST_P(CentralityGTest, testTopCloseness) {
    constexpr count size = 400;
    constexpr count k = 10;
    Aux::Random::setSeed(42, false);
    const auto G1 = DorogovtsevMendesGenerator(size).generate();
    Graph G(G1, false, isDirected());

    Closeness cc(G1, true, ClosenessVariant::generalized);
    cc.run();
    auto exactScores = cc.scores();
    auto ranking = cc.ranking();
    for (auto firstHeu : {true, false}) {
        for (auto secHeu : {true, false}) {
            TopCloseness topcc(G, k, firstHeu, secHeu);
            topcc.run();
            auto scores = topcc.topkScoresList();
            EXPECT_EQ(topcc.topkNodesList().size(), k);
            for (count i = 0; i < k; i++) {
                EXPECT_DOUBLE_EQ(ranking[i].second, scores[i]);
            }
        }
    }
}

TEST_P(CentralityGTest, testTopHarmonicCloseness) {
    const count size = 400;
    const double tol = 1e-6;

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(size, 0.01, isDirected()).generate();
        if (isWeighted()) {
            G = GraphTools::toWeighted(G);
            G.forEdges([&G](node u, node v) { G.setWeight(u, v, Aux::Random::probability()); });
        }
        HarmonicCloseness cc(G, false);
        cc.run();
        const auto ranking = cc.ranking();
        for (bool useNBbound : {true, false}) {
            for (count k : {5, 10, 20}) {
                if (isWeighted() && useNBbound)
                    continue;
                TopHarmonicCloseness topcc(G, k, useNBbound);
                topcc.run();

                auto topkScores = topcc.topkScoresList();
                EXPECT_EQ(topcc.topkNodesList().size(), k);
                EXPECT_EQ(topkScores.size(), k);

                topkScores = topcc.topkScoresList(true);

                for (count i = 0; i < topkScores.size(); ++i)
                    EXPECT_NEAR(ranking[i].second, topkScores[i], tol);
                for (count i = k; i < topkScores.size(); ++i)
                    EXPECT_NEAR(topkScores[i], topkScores[k - 1], tol);
            }
        }
    }
}

TEST_F(CentralityGTest, testLaplacianCentrality) {
    // The graph structure and reference values for the scores are taken from
    // Qi et al., Laplacian centrality: A new centrality measure for weighted
    // networks.
    //
    // See
    // https://math.wvu.edu/~cqzhang/Publication-files/my-paper/INS-2012-Laplacian-W.pdf.
    Graph G(6, true);

    G.addEdge(0, 1, 4);
    G.addEdge(0, 2, 2);
    G.addEdge(1, 2, 1);
    G.addEdge(1, 3, 2);
    G.addEdge(1, 4, 2);
    G.addEdge(4, 5, 1);

    LaplacianCentrality lc(G);
    lc.run();
    std::vector<double> scores = lc.scores();

    EXPECT_EQ(140, scores[0]);
    EXPECT_EQ(180, scores[1]);
    EXPECT_EQ(56, scores[2]);
    EXPECT_EQ(44, scores[3]);
    EXPECT_EQ(52, scores[4]);
    EXPECT_EQ(8, scores[5]);
}

TEST_F(CentralityGTest, testLaplacianCentralityNormalized) {
    Graph G(6, true);

    G.addEdge(0, 1, 4);
    G.addEdge(0, 2, 2);
    G.addEdge(1, 2, 1);
    G.addEdge(1, 3, 2);
    G.addEdge(1, 4, 2);
    G.addEdge(4, 5, 1);

    LaplacianCentrality lc(G, true);
    lc.run();
    std::vector<double> scores = lc.scores();

    EXPECT_EQ(0.70, scores[0]);
    EXPECT_EQ(0.90, scores[1]);
    EXPECT_EQ(0.28, scores[2]);
    EXPECT_EQ(0.22, scores[3]);
    EXPECT_EQ(0.26, scores[4]);
    EXPECT_EQ(0.04, scores[5]);
}

TEST_F(CentralityGTest, testLaplacianCentralityUnweighted) {
    Graph G(6);

    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(1, 4);
    G.addEdge(4, 5);

    LaplacianCentrality lc(G);
    lc.run();
    std::vector<double> scores = lc.scores();

    EXPECT_EQ(18, scores[0]);
    EXPECT_EQ(34, scores[1]);
    EXPECT_EQ(18, scores[2]);
    EXPECT_EQ(10, scores[3]);
    EXPECT_EQ(16, scores[4]);
    EXPECT_EQ(6, scores[5]);
}

TEST_P(CentralityGTest, testGroupDegree) {
    Aux::Random::setSeed(42, false);
    constexpr count nodes = 12;
    constexpr count k = 5;
    auto g = ErdosRenyiGenerator(nodes, 0.3, isDirected()).generate();

    auto computeGroupDegree = [&](const std::vector<bool> &curGroup, const Graph &g) {
        count result = 0;
        g.forNodes([&](node u) {
            if (!curGroup[u]) {
                bool neighborInGroup = false;
                g.forInNeighborsOf(u, [&](node v) {
                    if (!neighborInGroup && curGroup[v]) {
                        neighborInGroup = true;
                        ++result;
                    }
                });
            }
        });

        return result;
    };

    GroupDegree gd(g, k, false);
    gd.run();
    auto scoreNoGroup = gd.getScore();

    GroupDegree gdIncludeGroup(g, k, true);
    gdIncludeGroup.run();
    auto scorePlusGroup = gdIncludeGroup.getScore();

    std::vector<bool> reference(nodes, false);
    for (count i = nodes - k; i < nodes; ++i) {
        reference[i] = true;
    }

    count maxScore = 0;

    do {
        count curScore = computeGroupDegree(reference, g);
        if (curScore > maxScore) {
            maxScore = curScore;
        }
    } while (std::next_permutation(reference.begin(), reference.end()));

    EXPECT_TRUE(scoreNoGroup > 0.5 * maxScore);
    EXPECT_TRUE(scorePlusGroup >
                (1.0 - 1.0 / std::exp(1.0)) * static_cast<double>(maxScore + k));
    EXPECT_EQ(scoreNoGroup, gd.scoreOfGroup(gd.groupMaxDegree()));
    EXPECT_EQ(scorePlusGroup,
              gdIncludeGroup.scoreOfGroup(gdIncludeGroup.groupMaxDegree()));
}

TEST_F(CentralityGTest, testGroupBetweennessScore) {
/** 0           5
 *  | \       / |
 *  1 - 3 - 4 - 6
 *  | /       \ |
 *  2           7
 */
    Graph G(8, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(0, 3);
    G.addEdge(1, 3);
    G.addEdge(2, 3);

    G.addEdge(3, 4);

    G.addEdge(4, 5);
    G.addEdge(4, 6);
    G.addEdge(4, 7);
    G.addEdge(5, 6);
    G.addEdge(6, 7);

    BFS bfs(G, 0, true, false);
    // Naively computes the group betweenness of a group of nodes S
    auto naiveGB = [&](const std::vector<node> &S) -> double {
        double score = 0;
        count n = G.upperNodeIdBound();
        std::vector<bool> inGroup(n);
        for (node u : S)
            inGroup[u] = true;
        for (node source = 0; source < n; ++source) {
            bfs.setSource(source);
            bfs.run();
            for (node target = 0; target < n; ++target) {
                if (target == source)
                    continue;

                auto paths = bfs.getPaths(target);
                if (paths.empty())
                    continue;
                double curScore = 0;
                for (auto &path : paths) {
                    for (node u : path) {
                        if (u != source && u != target && inGroup[u]) {
                            curScore += 1;
                            break;
                        }
                    }
                }

                score += curScore / static_cast<double>(paths.size());
            }
        }

        return score;
    };

    ApproxGroupBetweenness gb(G, 2, 0.1);
    count z = G.numberOfNodes();
    for (count k = 1; k <= z; ++k) {
        std::vector<bool> inGroup(z);
        for (index i = 0; i < k; ++i)
            inGroup[z - i - 1] = true;
        do {
            std::vector<node> group;
            group.reserve(k);
            for (node u = 0; u < z; ++u) {
                if (inGroup[u]) {
                    group.push_back(u);
                    if (group.size() == k)
                        break;
                }
            }
            double computedScore = gb.scoreOfGroup(group);
            double naiveScore = naiveGB(group);
            EXPECT_NEAR(computedScore, naiveScore, 1e-6);
        } while (std::next_permutation(inGroup.begin(), inGroup.end()));
    }
}

TEST_F(CentralityGTest, runTestApproxGroupBetweennessSmallGraph) {

    Aux::Random::setSeed(42, false);

    count n = 8;
    Graph g(n, false, false);

    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 4);
    g.addEdge(3, 5);
    g.addEdge(4, 5);
    g.addEdge(5, 6);
    g.addEdge(5, 7);
    g.addEdge(0, 5);

    count k = 2;
    double eps = 0.1;
    ApproxGroupBetweenness gb(g, k, eps);
    gb.run();

    double maxScore = 0;
    std::vector<bool> inGroup(n);
    for (index i = 0; i < k; ++i)
        inGroup[n - i - 1] = true;
    do {
        std::vector<node> group;
        group.reserve(k);
        for (node u = 0; u < n; ++u) {
            if (inGroup[u]) {
                group.push_back(u);
                if (group.size() == k)
                    break;
            }
        }

        maxScore = std::max(maxScore, gb.scoreOfGroup(group));
    } while (std::next_permutation(inGroup.begin(), inGroup.end()));

    EXPECT_TRUE(gb.scoreOfGroup(gb.groupMaxBetweenness()) >= maxScore * eps);
}

TEST_F(CentralityGTest, testGroupCloseness) {
    Aux::Random::setSeed(42, false);

    Graph g(8, false, false);

    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 4);
    g.addEdge(3, 5);
    g.addEdge(4, 5);
    g.addEdge(5, 6);
    g.addEdge(5, 7);
    g.addEdge(0, 5);

    count k = 3;

    GroupCloseness gc(g, k);
    gc.run();
    auto apx = gc.groupMaxCloseness();
    std::sort(apx.begin(), apx.end());
    std::vector<node> solution = {0, 2, 5};
    for (count i = 0; i < k; ++i) {
        EXPECT_EQ(apx[i], solution[i]);
    }

    EXPECT_NEAR(gc.scoreOfGroup(apx), 1.0, 1e-5);
}

/**
 * This test succeeds with the fixed random seed (42).
 * However, the Kadabra algorithm computes a correct epsilon-approximation of
 * the betweenness centrality score of all the nodes of the graph with high
 * probability. Thus, it is possible that, for a different random seed, this
 * test fails.
 */
TEST_F(CentralityGTest, testKadabraAbsolute) {
    Aux::Random::setSeed(42, true);
    const count n = 10;
    Graph g = ErdosRenyiGenerator(n, 0.1).generate();

    const double delta = 0.1;
    const double epsilon = 0.01;
    KadabraBetweenness kadabra(g, epsilon, delta);
    kadabra.run();
    auto scores = kadabra.topkScoresList();
    auto nodes = kadabra.topkNodesList();

    Betweenness betweenness(g, true);
    betweenness.run();
    count maxErrors = (count)std::ceil(delta * (double)n);

    count errors = 0;
    for (count i = 0; i < n; ++i) {
        if (std::abs(scores[i] - betweenness.score(nodes[i])) > delta) {
            ++errors;
        }
    }

    EXPECT_TRUE(errors <= maxErrors);
}

/**
 * This test succeeds with the fixed random seed (42).
 * However, the Kadabra algorithm finds the top-k nodes with
 * highest betweenness centrality with high probability. Thus, it is possible
 * that, for a different random seed, this test fails.
 */

TEST_F(CentralityGTest, testKadabraTopK) {
    Aux::Random::setSeed(42, true);
    const count n = 10;
    Graph g = ErdosRenyiGenerator(n, 0.1).generate();

    const double delta = 0.1;
    const double epsilon = 0.01;
    const count k = 3;
    KadabraBetweenness kadabra(g, epsilon, delta, k);
    kadabra.run();
    auto kadabraRanking = kadabra.ranking();

    Betweenness betweenness(g, true);
    betweenness.run();
    auto betwRanking = betweenness.ranking();
    bool correctRanking = true;
    for (count i = 0; i < k; ++i) {
        if (betwRanking[i].first != kadabraRanking[i].first) {
            correctRanking = false;
            int j = static_cast<int>(i) - 1;
            while (j >= 0 && betwRanking[j].second == betwRanking[i].second) {
                --j;
            }
            ++j;
            while (j < static_cast<int>(n) && betwRanking[j].second == betwRanking[i].second) {
                if (betwRanking[j].first == kadabraRanking[i].first) {
                    correctRanking = true;
                    break;
                }
                ++j;
            }
        }
    }
    EXPECT_TRUE(correctRanking);
}

TEST_F(CentralityGTest, testKadabraAbsoluteDeterministic) {
    const index seed = 42;
    Aux::Random::setSeed(seed, true);
    const count n = 300;
    Graph g = ErdosRenyiGenerator(n, 0.1).generate();

    const double delta = 0.1;
    const double epsilon = 0.01;
    Aux::Random::setSeed(seed, true);
    KadabraBetweenness kadabra(g, epsilon, delta, true);
    kadabra.run();
    auto scores = kadabra.topkScoresList();
    auto nodes = kadabra.topkNodesList();

    Aux::Random::setSeed(seed, true);
    KadabraBetweenness kadabra2(g, epsilon, delta, true);
    kadabra2.run();
    auto scores2 = kadabra2.topkScoresList();
    auto nodes2 = kadabra2.topkNodesList();

    for (count i = 0; i < n; ++i) {
        EXPECT_EQ(scores[i], scores2[i]);
    }
}

TEST_P(CentralityGTest, testDynTopHarmonicCloseness) {
    auto G1 = DorogovtsevMendesGenerator(500).generate();
    Graph G(G1, false, isDirected());

    constexpr count k = 10;

    DynTopHarmonicCloseness centrality(G, k, false);
    centrality.run();

    HarmonicCloseness reference(G, false);
    reference.run();

    auto scores = centrality.ranking();
    auto refScores = reference.ranking();

    for (count j = 0; j < k; ++j) {
        EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
    }

    count numInsertions = 1;

    std::vector<GraphEvent> deletions;
    std::vector<GraphEvent> insertions;

    for (count i = 0; i < numInsertions; i++) {

        node u = G.upperNodeIdBound();
        node v = G.upperNodeIdBound();

        do {
            u = GraphTools::randomNode(G);
            v = GraphTools::randomNode(G);
        } while (G.hasEdge(u, v));

        GraphEvent edgeAddition(GraphEvent::EDGE_ADDITION, u, v);
        insertions.insert(insertions.begin(), edgeAddition);

        GraphEvent edgeDeletion(GraphEvent::EDGE_REMOVAL, u, v);
        deletions.push_back(edgeDeletion);

        G.addEdge(u, v);
    }

    for (auto &e : insertions) {
        G.removeEdge(e.u, e.v);
    }

    for (auto &edgeAddition : insertions) {

        node u = edgeAddition.u;
        node v = edgeAddition.v;

        G.addEdge(u, v);

        centrality.update(edgeAddition);
        reference.run();

        scores = centrality.ranking();
        refScores = reference.ranking();

        for (count j = 0; j < k; ++j) {
            EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
        }
    }

    for (auto &edgeDeletion : deletions) {
        node u = edgeDeletion.u;
        node v = edgeDeletion.v;

        G.removeEdge(u, v);

        centrality.update(edgeDeletion);
        reference.run();

        scores = centrality.ranking();
        refScores = reference.ranking();

        for (count j = 0; j < k; ++j) {
            EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
        }
    }

    for (auto &edgeInsertion : insertions) {
        G.addEdge(edgeInsertion.u, edgeInsertion.v);
    }

    reference.run();
    centrality.updateBatch(insertions);

    scores = centrality.ranking();
    refScores = reference.ranking();

    for (count j = 0; j < k; ++j) {
        EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
    }
}

TEST_F(CentralityGTest, testApproxSpanningEdge) {
    Aux::Random::setSeed(42, false);

    // Testing a graph that is too small yields approximation errors; however, if the input
    // graph has > 100 nodes, running the exact algorithm takes too long.
    // Therefore, in this test we use the approximation algorithm implemented in
    // SpanningEdgeCentrality as baseline, and check that the values are within a 2-epsilon
    // approximation.
    Graph G = ErdosRenyiGenerator(300, 0.1, false).generate();
    G.indexEdges();
    constexpr double eps = 0.1;

    ApproxSpanningEdge apx(G, eps);
    apx.run();
    SpanningEdgeCentrality se(G, eps);
    se.runParallelApproximation();
    auto apxScores = apx.scores();
    auto exactScores = se.scores();

    G.forEdges([&](node /*u*/, node /*v*/, edgeweight /*w*/, edgeid eid) {
        EXPECT_NEAR(apxScores[eid], exactScores[eid], 2*eps);
    });
}

TEST_P(CentralityGTest, testGedWalk) {
    Aux::Random::setSeed(42, true);
    constexpr count k = 3;
    constexpr double epsilon = 0.01;
    auto g = ErdosRenyiGenerator(20, 0.1, isDirected()).generate();

    for (const auto bs : {GedWalk::BoundStrategy::geometric, GedWalk::BoundStrategy::spectral}) {
        for (const auto gs : {GedWalk::GreedyStrategy::lazy, GedWalk::GreedyStrategy::stochastic}) {
            GedWalk gedWalk(g, k, epsilon, -1.0, bs, gs);
            gedWalk.run();

            const auto apxScore = gedWalk.getApproximateScore();
            EXPECT_GE(apxScore, 0);
            const auto apxGroup = gedWalk.groupMaxGedWalk();
            EXPECT_EQ(std::unordered_set<node>(apxGroup.begin(), apxGroup.end()).size(), k);

            double maxScore = 0.;
            std::vector<node> group;
            std::vector<node> max_group;
            std::vector<bool> reference(g.numberOfNodes(), false);
            std::fill(reference.end() - k, reference.end(), true);

            do {
                group.clear();
                for (count i = 0; i < reference.size(); ++i) {
                    if (reference[i]) {
                        group.push_back(i);
                    }
                }

                const auto curScore = gedWalk.scoreOfGroup(group.begin(), group.end(), epsilon);
                if (curScore > maxScore) {
                    maxScore = curScore;
                    max_group = group;
                }
            } while (std::next_permutation(reference.begin(), reference.end()));

            EXPECT_GE(apxScore, (1. - 1. / std::exp(1)) * maxScore - epsilon);
        }
    }
}

TEST_F(CentralityGTest, testApproxElectricalCloseness) {
    const double eps = 0.1;
    const count n = 75;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);
        auto G = HyperbolicGenerator(n, 10, 3).generate();
        G = ConnectedComponents::extractLargestConnectedComponent(G, true);

        // Create a biconnected component with size 2.
        G.addNodes(2);
        G.addEdge(n - 1, n);
        G.addEdge(n, n + 1);

        ApproxElectricalCloseness apx(G);
        apx.run();
        const auto diag = apx.getDiagonal();
        const auto gt = apx.computeExactDiagonal(1e-12);
        G.forNodes([&](node u) { EXPECT_NEAR(diag[u], gt[u], eps); });
        EXPECT_EQ(apx.scores().size(), G.numberOfNodes());
    }
}

TEST_P(CentralityGTest, testGroupClosenessGrowShrink) {
    if (isDirected()) { // directed graphs are not supported
        Graph G(10, isWeighted(), true);
        std::array<node, 1> group;
        EXPECT_THROW(GroupClosenessGrowShrink(G, group.begin(), group.end()), std::runtime_error);
        return;
    }

    { // Empty input groups are not supported
        Graph G(10, isWeighted(), false);
        std::vector<node> emptyGroup;
        EXPECT_THROW(GroupClosenessGrowShrink(G, emptyGroup.begin(), emptyGroup.end()),
                     std::runtime_error);
    }

    const count k = 5;
    auto G = EdgeListReader{'\t', 0, "#", false, false}.read("input/MIT8.edgelist");
    G = ConnectedComponents::extractLargestConnectedComponent(G);

    if (isWeighted()) {
        G = GraphTools::toWeighted(G);
        G.forEdges([&G](const node u, const node v) {
            G.setWeight(u, v, Aux::Random::probability());
        });
    }

    auto farnessOfGroup = [&](const Graph &G, const std::unordered_set<node> &group) -> edgeweight {
        edgeweight farness = 0;
        if (G.isWeighted()) {
            Traversal::DijkstraFrom(
                G, group.begin(), group.end(),
                [&farness](node, const edgeweight distance) { farness += distance; });
        } else {
            Traversal::BFSfrom(
                G, group.begin(), group.end(),
                [&farness](node, const edgeweight distance) { farness += distance; });
        }

        return farness;
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);

        std::unordered_set<node> group;

        do {
            group.insert(GraphTools::randomNode(G));
        } while (group.size() < k);

        edgeweight sumDist = farnessOfGroup(G, group);
        GroupClosenessGrowShrink gc(G, group.begin(), group.end());
        gc.run();
        auto groupMaxCC = gc.groupMaxCloseness();
        count nSwaps = gc.numberOfIterations();

        EXPECT_EQ(groupMaxCC.size(), k);

        edgeweight sumDistGroupMaxCC =
            farnessOfGroup(G, std::unordered_set<node>(groupMaxCC.begin(), groupMaxCC.end()));
        if (nSwaps > 0) {
            EXPECT_GE(sumDist, sumDistGroupMaxCC);
        } else {
            EXPECT_EQ(sumDist, sumDistGroupMaxCC);
            std::for_each(groupMaxCC.begin(), groupMaxCC.end(), [&group](const node u) {
                EXPECT_NE(group.find(u), group.end());
            });
        }
    }
}


TEST_P(CentralityGTest, testDegreeCentrality) {
    Graph g(8, false, isDirected());

    g.addEdge(0, 2);
    g.addEdge(0, 5);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 2);
    g.addEdge(2, 4);
    g.addEdge(3, 5);
    g.addEdge(4, 5);
    g.addEdge(5, 5);
    g.addEdge(5, 6);
    g.addEdge(5, 7);
    g.addEdge(7, 7);

    DegreeCentrality dc(g, false, true, false);
    dc.run();

    if (isDirected()){
        std::array<int ,8> expectedResults{{2, 1, 3, 1, 1, 3, 0, 1}};
        for(long unsigned int i = 0; i < expectedResults.size(); i++) EXPECT_EQ(expectedResults[i], dc.score(i));
    } else {
        std::array<int ,8> expectedResults{{2, 1, 5, 2, 2, 6, 1, 2}};
        for(long unsigned int i = 0; i < expectedResults.size(); i++) EXPECT_EQ(expectedResults[i], dc.score(i));
    }
}

TEST_P(CentralityGTest, testInDegreeCentrality) {
    Graph g(8, false, isDirected());

    g.addEdge(0, 2);
    g.addEdge(0, 5);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 2);
    g.addEdge(2, 4);
    g.addEdge(3, 5);
    g.addEdge(4, 5);
    g.addEdge(5, 5);
    g.addEdge(5, 6);
    g.addEdge(5, 7);
    g.addEdge(7, 7);

    DegreeCentrality dc(g, false, false, false);
    dc.run();

    if (isDirected()){
        std::array<int ,8> expectedResults{{0, 0, 3, 1, 1, 4, 1, 2}};
        for(long unsigned int i = 0; i < expectedResults.size(); i++) EXPECT_EQ(expectedResults[i], dc.score(i));
    } else {
        std::array<int ,8> expectedResults{{2, 1, 5, 2, 2, 6, 1, 2}};
        for(long unsigned int i = 0; i < expectedResults.size(); i++) EXPECT_EQ(expectedResults[i], dc.score(i));
    }
}

TEST_P(CentralityGTest, testDegreeCentralityIgnoreSelfLoops) {
    Graph g(8, false, isDirected());

    g.addEdge(0, 2);
    g.addEdge(0, 5);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 2);
    g.addEdge(2, 4);
    g.addEdge(3, 5);
    g.addEdge(4, 5);
    g.addEdge(5, 5);
    g.addEdge(5, 6);
    g.addEdge(5, 7);
    g.addEdge(7, 7);

    DegreeCentrality dc(g, false, true, true);
    dc.run();

    if (isDirected()){
        std::array<int ,8> expectedResults{{2, 1, 2, 1, 1, 2, 0, 0}};
        for(long unsigned int i = 0; i < expectedResults.size(); i++) EXPECT_EQ(expectedResults[i], dc.score(i));
    } else {
        std::array<int ,8> expectedResults{{2, 1, 4, 2, 2, 5, 1, 1}};
        for(long unsigned int i = 0; i < expectedResults.size(); i++) EXPECT_EQ(expectedResults[i], dc.score(i));
    }
}

TEST_P(CentralityGTest, testGroupClosenessLocalSwaps) {
    if (isDirected()) { // directed graphs are not supported
        Graph G(10, isWeighted(), true);
        std::array<node, 1> group;
        EXPECT_THROW(GroupClosenessLocalSwaps(G, group.begin(), group.end()), std::runtime_error);
        return;
    }

    { // Empty input groups are not supported
        Graph G(10, isWeighted(), false);
        std::vector<node> emptyGroup;
        EXPECT_THROW(GroupClosenessLocalSwaps(G, emptyGroup.begin(), emptyGroup.end()), std::runtime_error);
    }

    const count k = 5;
    auto G = EdgeListReader{'\t', 0, "#", false, false}.read("input/MIT8.edgelist");
    G = ConnectedComponents::extractLargestConnectedComponent(G);

    if (isWeighted()) {
        G = GraphTools::toWeighted(G);
        G.forEdges([&G](const node u, const node v) {
            G.setWeight(u, v, Aux::Random::probability());
        });
    }

    auto farnessOfGroup = [&](const Graph &G, const std::unordered_set<node> &group) -> count {
        count farness = 0;
        Traversal::BFSfrom(G, group.begin(), group.end(),
                           [&farness](node, count distance) { farness += distance; });

        return farness;
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);

        std::unordered_set<node> group;

        do {
            group.insert(GraphTools::randomNode(G));
        } while (group.size() < k);

        const count maxSwaps = 100;
        const count sumDist = farnessOfGroup(G, group);
        GroupClosenessLocalSwaps gc(G, group.begin(), group.end(), maxSwaps);
        gc.run();
        auto groupMaxCC = gc.groupMaxCloseness();
        const count nSwaps = gc.numberOfSwaps();

        EXPECT_LE(nSwaps, maxSwaps);
        EXPECT_EQ(groupMaxCC.size(), k);

        count sumDistGroupMaxCC =
            farnessOfGroup(G, std::unordered_set<node>(groupMaxCC.begin(), groupMaxCC.end()));

        if (nSwaps > 0) {
            EXPECT_GE(sumDist, sumDistGroupMaxCC);
        } else {
            EXPECT_EQ(sumDist, sumDistGroupMaxCC);
            std::for_each(groupMaxCC.begin(), groupMaxCC.end(), [&group](const node u) {
                EXPECT_NE(group.find(u), group.end());
            });

        }
    }
}

TEST_P(CentralityGTest, testGroupHarmonicCloseness) {

    const auto computeOpt = [&](const Graph &G, count k) -> double {
        std::vector<bool> inGroup(G.upperNodeIdBound());
        std::fill(inGroup.begin(), inGroup.begin() + k, true);
        double opt = -std::numeric_limits<double>::max();
        std::vector<node> group;
        group.reserve(k);
        do {
            group.clear();
            G.forNodes([&](node u) {
                if (inGroup[u])
                    group.push_back(u);
            });
            opt = std::max(opt, GroupHarmonicCloseness::scoreOfGroup(G, group.begin(), group.end()));
        } while (std::prev_permutation(inGroup.begin(), inGroup.end()));

        return opt;
    };

    const double weightUB = 10;
    for (const int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);
        auto G = ErdosRenyiGenerator(20, 0.1, isDirected()).generate();
        double lambda = 1;

        if (isWeighted()) {
            double maxWeight = 0, minWeight = weightUB;
            G = GraphTools::toWeighted(G);
            G.forEdges([&](node u, node v) {
                const double curWeight = Aux::Random::real(0.001, weightUB);
                maxWeight = std::max(maxWeight, curWeight);
                minWeight = std::min(minWeight, curWeight);
                G.setWeight(u, v, curWeight);
            });
            lambda = minWeight / maxWeight;
        }

        // Guaranted approximation ratio
        double approxRatio;
        if (isDirected())
            approxRatio = lambda * (1. - 1. / (2. * std::exp(1.)));
        else
            approxRatio = lambda * (1. - 1. / std::exp(1.)) / 2.;

        for (const count k : {3, 4, 5}) {
            GroupHarmonicCloseness ghc(G, k);
            ghc.run();
            const auto group = ghc.groupMaxHarmonicCloseness();

            // Test group
            EXPECT_EQ(group.size(), k);
            EXPECT_EQ(std::unordered_set<node>(group.begin(), group.end()).size(), k);

            // Test quality
            const double score = GroupHarmonicCloseness::scoreOfGroup(G, group.begin(), group.end());
            const double opt = computeOpt(G, k);
            EXPECT_GE(opt, score);
            EXPECT_GE(score / opt, approxRatio);
        }
    }
}

TEST_P(CentralityGTest, testGroupClosenessLocalSearch) {
    { // Empty groups are not allowed
        std::vector<node> emptyVector;
        Graph G(10, isWeighted(), isDirected());
        EXPECT_THROW(GroupClosenessLocalSearch(G, emptyVector.begin(), emptyVector.end()),
                     std::runtime_error);
    }

    const auto groupCloseness = [&](const Graph &G, const std::vector<node> &group) -> edgeweight {
        edgeweight groupFarness = 0;
        Traversal::DijkstraFrom(G, group.begin(), group.end(),
                                [&groupFarness](node, edgeweight dist) { groupFarness += dist; });
        return groupFarness > 0 ? 1. / groupFarness : 0;
    };

    Aux::Random::setSeed(1, true);
    auto G = ErdosRenyiGenerator(100, 0.1, isDirected()).generate();
    if (isWeighted()) {
        G = GraphTools::toWeighted(G);
        G.forEdges([&](node u, node v) { G.setWeight(u, v, Aux::Random::real(10)); });
    }

    std::unordered_set<node> initGroup;
    const count k = 5;
    do {
        initGroup.insert(GraphTools::randomNode(G));
    } while (initGroup.size() < k);

    const auto initGC = groupCloseness(G, std::vector<node>{initGroup.begin(), initGroup.end()});

    GroupClosenessLocalSearch gcls(G, initGroup.begin(), initGroup.end(), !G.isDirected());
    gcls.run();

    const auto group = gcls.groupMaxCloseness();
    EXPECT_EQ(std::unordered_set<node>(group.begin(), group.end()).size(), k);
    EXPECT_GE(groupCloseness(G, group), initGC);

    GroupClosenessLocalSearch gcls2(G, group, false);
    gcls2.run();

    EXPECT_EQ(gcls2.numberOfIterations(), 0);
}

} /* namespace NetworKit */
