// no-networkit-format
/*
 * AlgebraicTriangleCountingGTest.cpp
 *
 *  Created on: Jul 12, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <gtest/gtest.h>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/algorithms/AlgebraicTriangleCounting.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/centrality/LocalClusteringCoefficient.hpp>

namespace NetworKit {

class AlgebraicTriangleCountingGTest : public testing::Test {};

TEST(AlgebraicTriangleCountingGTest, testToyGraphOne) {
    Graph graph(5);

    graph.addEdge(0,1);
    graph.addEdge(0,2);
    graph.addEdge(1,2);

    AlgebraicTriangleCounting<CSRMatrix> atc(graph);
    atc.run();

    std::vector<count> nodeScores = atc.getScores();

    EXPECT_EQ(1u, nodeScores[0]) << "wrong triangle count";
    EXPECT_EQ(1u, nodeScores[1]) << "wrong triangle count";
    EXPECT_EQ(1u, nodeScores[2]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeScores[3]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeScores[4]) << "wrong triangle count";
}

TEST(AlgebraicTriangleCountingGTest, testToyGraphTwo) {
    Graph graph(5);

    graph.addEdge(0,1);
    graph.addEdge(0,2);
    graph.addEdge(1,2);
    graph.addEdge(1,3);
    graph.addEdge(3,4);
    graph.addEdge(2,4);

    AlgebraicTriangleCounting<CSRMatrix> atc(graph);
    atc.run();

    EXPECT_TRUE(atc.hasFinished());
    std::vector<count> nodeScores = atc.getScores();

    EXPECT_EQ(1u, nodeScores[0]) << "wrong triangle count";
    EXPECT_EQ(1u, nodeScores[1]) << "wrong triangle count";
    EXPECT_EQ(1u, nodeScores[2]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeScores[3]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeScores[4]) << "wrong triangle count";

    graph.addEdge(2,3);
    atc = AlgebraicTriangleCounting<CSRMatrix>(graph);
    atc.run();

    EXPECT_EQ(1u, atc.score(0)) << "wrong triangle count";
    EXPECT_EQ(2u, atc.score(1)) << "wrong triangle count";
    EXPECT_EQ(3u, atc.score(2)) << "wrong triangle count";
    EXPECT_EQ(2u, atc.score(3)) << "wrong triangle count";
    EXPECT_EQ(1u, atc.score(4)) << "wrong triangle count";
}

TEST(AlgebraicTriangleCountingGTest, testDirectedToyGraphOne) {
    Graph graph(5, false, true);

    graph.addEdge(0,1);
    graph.addEdge(1,2);
    graph.addEdge(2,0);

    AlgebraicTriangleCounting<CSRMatrix> atc(graph);
    atc.run();

    std::vector<count> nodeScores = atc.getScores();

    EXPECT_EQ(1u, nodeScores[0]) << "wrong triangle count";
    EXPECT_EQ(1u, nodeScores[1]) << "wrong triangle count";
    EXPECT_EQ(1u, nodeScores[2]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeScores[3]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeScores[4]) << "wrong triangle count";
}

TEST(AlgebraicTriangleCountingGTest, testDirectedToyGraphTwo) {
    Graph graph(5, false, true);

    graph.addEdge(0,1);
    graph.addEdge(0,2);
    graph.addEdge(1,2);

    AlgebraicTriangleCounting<CSRMatrix> atc(graph);
    atc.run();

    std::vector<count> nodeScores = atc.getScores();

    EXPECT_EQ(0u, nodeScores[0]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeScores[1]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeScores[2]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeScores[3]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeScores[4]) << "wrong triangle count";
}

TEST(AlgebraicTriangleCountingGTest, testDirectedToyGraphThree) {
    Graph graph(5, false, true);

    graph.addEdge(0,1);
    graph.addEdge(1,0);
    graph.addEdge(0,2);
    graph.addEdge(2,0);
    graph.addEdge(1,2);
    graph.addEdge(2,1);

    AlgebraicTriangleCounting<CSRMatrix> atc(graph);
    atc.run();

    std::vector<count> nodeScores = atc.getScores();

    EXPECT_EQ(2u, nodeScores[0]) << "wrong triangle count";
    EXPECT_EQ(2u, nodeScores[1]) << "wrong triangle count";
    EXPECT_EQ(2u, nodeScores[2]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeScores[3]) << "wrong triangle count";
    EXPECT_EQ(0u, nodeScores[4]) << "wrong triangle count";
}

TEST(AlgebraicTriangleCountingGTest, testLocalClusteringCoefficient) {
    METISGraphReader reader;
    Graph graph = reader.read("input/celegans_metabolic.graph");
    INFO("graph has ", graph.numberOfNodes(), " nodes and ", graph.numberOfEdges(), " edges and directed? ", graph.isDirected());

    Aux::Timer timer;
    AlgebraicTriangleCounting<CSRMatrix> atc(graph);
    timer.start();
    atc.run();

    std::vector<count> nodeScores = atc.getScores();
    std::vector<double> lccAlgebraic(graph.numberOfNodes());
    graph.parallelForNodes([&](node u) {
        if (graph.degree(u) < 2) {
            lccAlgebraic[u] = 0.0;
        } else {
            lccAlgebraic[u] = 2.0 * nodeScores[u] / (graph.degree(u) * (graph.degree(u)-1));
        }
    });
    timer.stop();

    INFO("Algebraic local clustering coefficient took ", timer.elapsedTag());

    LocalClusteringCoefficient lcc(graph);
    timer.start();
    lcc.run();
    std::vector<double> lccValues = lcc.scores();
    timer.stop();

    INFO("Graph theoretic local clustering coefficient took ", timer.elapsedTag());

    graph.forNodes([&](node u) {
        EXPECT_EQ(lccValues[u], lccAlgebraic[u]);
    });
}

} /* namespace NetworKit */
