// no-networkit-format
/*
 * AlgebraicSpanningEdgeCentralityGTest.cpp
 *
 *  Created on: Jul 12, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */


#include <gtest/gtest.h>

#include <networkit/algebraic/algorithms/AlgebraicSpanningEdgeCentrality.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>

#include <networkit/io/METISGraphReader.hpp>
#include <networkit/centrality/SpanningEdgeCentrality.hpp>

namespace NetworKit {

class AlgebraicSpanningEdgeCentralityGTest : public testing::Test{};


TEST_F(AlgebraicSpanningEdgeCentralityGTest, testOnToyGraph) {
    /* Graph:
            0    3
             \  / \
              2    5
             /  \ /
            1    4
     */
    count n = 6;
    Graph G(n, false, false);
    G.indexEdges();


    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    AlgebraicSpanningEdgeCentrality<CSRMatrix> asp(G);

    asp.run();
    EXPECT_NEAR(1.0, asp.score(0), 1e-5);
    EXPECT_NEAR(1.0, asp.score(1), 1e-5);
    EXPECT_NEAR(0.75, asp.score(2), 1e-5);
    EXPECT_NEAR(0.75, asp.score(3), 1e-5);
    EXPECT_NEAR(0.75, asp.score(4), 1e-5);
    EXPECT_NEAR(0.75, asp.score(5), 1e-5);
}

TEST_F(AlgebraicSpanningEdgeCentralityGTest, benchmarkSpanning) {
    METISGraphReader reader;

    std::string graphFiles[2] = {"input/karate.graph", "input/tiny_01.graph"};

    for (auto graphFile: graphFiles) {
        Graph G = reader.read(graphFile);
        G.indexEdges();
        Aux::Timer timer;

        SpanningEdgeCentrality exact(G);
        timer.start();
        exact.run();
        timer.stop();
        INFO("exact spanning edge centrality ranking time: ", timer.elapsedTag());

        AlgebraicSpanningEdgeCentrality<CSRMatrix> asp(G);

        timer.start();
        asp.runApproximation();
        timer.stop();
        INFO("algebraic spanning edge centrality ranking time: ", timer.elapsedTag());

        double error = 0.0;
        G.forEdges([&](node, node, edgeid e) {
            double relError = fabs(asp.score(e) - exact.score(e));
            if (fabs(exact.score(e)) > 1e-9) relError /= exact.score(e);
            error += relError;
        });
        error /= G.numberOfEdges();
        INFO("Avg. relative error: ", error);

        SpanningEdgeCentrality sp(G);
        timer.start();
        sp.runParallelApproximation();
        timer.stop();
        INFO("graph theoretical spanning edge centrality time: ", timer.elapsedTag());



        error = 0.0;
        G.forEdges([&](node, node, edgeid e) {
            double relError = fabs(sp.score(e) - exact.score(e));
            if (fabs(exact.score(e)) > 1e-9) relError /= exact.score(e);
            error += relError;
        });
        error /= G.numberOfEdges();
        INFO("Avg. relative error: ", error);
    }
}

} /* namespace NetworKit */
