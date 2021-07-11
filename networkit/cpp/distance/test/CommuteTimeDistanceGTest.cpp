// no-networkit-format
/*
 * CommuteTimeDistanceGTest.cpp
 *
 *  Created on: Jan 17, 2016
 *      Author: Michael
 */

#include <gtest/gtest.h>

#include <networkit/distance/CommuteTimeDistance.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/centrality/SpanningEdgeCentrality.hpp>

#include <tlx/unused.hpp>

#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>


namespace NetworKit {
class CommuteTimeDistanceGTest : public testing::Test{};

TEST_F(CommuteTimeDistanceGTest, testOnToyGraph) {
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

    SpanningEdgeCentrality sp(G);

    sp.run();
    EXPECT_NEAR(1.0, sp.score(0), 1e-5);
    EXPECT_NEAR(1.0, sp.score(1), 1e-5);
    EXPECT_NEAR(0.75, sp.score(2), 1e-5);
    EXPECT_NEAR(0.75, sp.score(3), 1e-5);
    EXPECT_NEAR(0.75, sp.score(4), 1e-5);
    EXPECT_NEAR(0.75, sp.score(5), 1e-5);

    CommuteTimeDistance ctd(G);
    ctd.run();
    double volG = 2.0 * G.numberOfEdges();
    EXPECT_NEAR(sqrt(1.0 * volG), ctd.distance(0, 2), 1e-4);
    EXPECT_NEAR(sqrt(1.0 * volG), ctd.distance(1, 2), 1e-4);
    EXPECT_NEAR(sqrt(0.75 * volG), ctd.distance(2, 3), 1e-4);
    EXPECT_NEAR(sqrt(0.75 * volG), ctd.distance(2, 4), 1e-4);
    EXPECT_NEAR(sqrt(0.75 * volG), ctd.distance(3, 5), 1e-4);
    EXPECT_NEAR(sqrt(0.75 * volG), ctd.distance(4, 5), 1e-4);
}

TEST_F(CommuteTimeDistanceGTest, testOnWeightedToyGraph) {
    /* Graph:
            0    3
             \  / \
              2    5
             /  \ /
            1    4

                Edge weight of (u,v): u+v
     */
    count n = 6;
    Graph G(n, true, false);
    G.indexEdges();


    G.addEdge(0, 2, 2);
    G.addEdge(1, 2, 3);
    G.addEdge(2, 3, 5);
    G.addEdge(2, 4, 6);
    G.addEdge(3, 5, 8);
    G.addEdge(4, 5, 9);

    CommuteTimeDistance ctd(G);
    ctd.run();
    double volG = 2.0 * G.totalEdgeWeight();
    DEBUG("volume : ", volG);
    EXPECT_NEAR(sqrt(0.5 * volG), ctd.distance(0, 2), 1e-3);
    EXPECT_NEAR(sqrt(0.3333 * volG), ctd.distance(1, 2), 1e-3);
    EXPECT_NEAR(sqrt(0.1336 * volG), ctd.distance(2, 3), 1e-3);
    EXPECT_NEAR(sqrt(0.1206 * volG), ctd.distance(2, 4), 1e-3);
    EXPECT_NEAR(sqrt(0.0991 * volG), ctd.distance(3, 5), 1e-3);
    EXPECT_NEAR(sqrt(0.0906 * volG), ctd.distance(4, 5), 1e-3);
}

TEST_F(CommuteTimeDistanceGTest, runECTDOnSmallGraphs) {
    METISGraphReader reader;

    std::string graphFiles[1] = { "input/tiny_01.graph"};

    for (auto graphFile: graphFiles) {
        Graph G = reader.read(graphFile);
        G.indexEdges();
        Aux::Timer timer;
        CommuteTimeDistance exact(G);
        CommuteTimeDistance cen(G);

        timer.start();
        exact.run();
        timer.stop();
        DEBUG("ECTD time: ", timer.elapsedTag());

        timer.start();
        cen.runApproximation();
        timer.stop();
        DEBUG("approx ECTD time: ", timer.elapsedTag());

        double error = 0.0;
        G.forNodes([&](node u){
            G.forNodes([&](node v) {
                double relError = fabs(cen.distance(u,v) - exact.distance(u,v));
            //	INFO("Approximated: ", cen.distance(u,v), ", exact: ", exact.distance(u,v));
                if (fabs(exact.distance(u,v)) > 1e-9) {
                    relError /= exact.distance(u,v);
                }
                error += relError;
            });
        });
        error /= G.numberOfNodes()*G.numberOfNodes();
        DEBUG("Avg. relative error: ", error);
    }
}

TEST_F(CommuteTimeDistanceGTest, runECTDParallelOnSmallGraphs) {
    METISGraphReader reader;

    std::string graphFiles[1] = { "input/tiny_01.graph"};

    for (auto graphFile: graphFiles) {
        Graph G = reader.read(graphFile);
        G.indexEdges();
        Aux::Timer timer;
        CommuteTimeDistance exact(G);
        CommuteTimeDistance cen(G);

        timer.start();
        exact.run();
        timer.stop();
        DEBUG("ECTD time: ", timer.elapsedTag());

        timer.start();
        cen.runParallelApproximation();
        timer.stop();
        DEBUG("approx ECTD time: ", timer.elapsedTag());

        double error = 0.0;
        G.forNodes([&](node u){
            G.forNodes([&](node v) {
                double relError = fabs(cen.distance(u,v) - exact.distance(u,v));
            //	INFO("Approximated: ", cen.distance(u,v), ", exact: ", exact.distance(u,v));
                if (fabs(exact.distance(u,v)) > 1e-9) {
                    relError /= exact.distance(u,v);
                }
                error += relError;
            });
        });
        error /= G.numberOfNodes()*G.numberOfNodes();
        DEBUG("Avg. relative error: ", error);
    }
}

TEST_F(CommuteTimeDistanceGTest, runECTDSingleSource) {
    METISGraphReader reader;

    std::string graphFiles[2] = {"input/karate.graph", "input/tiny_01.graph"};

    for (auto graphFile: graphFiles) {
        Graph G = reader.read(graphFile);
        CommuteTimeDistance ectd(G);
        node u = GraphTools::randomNode(G);
        double sum1 = ectd.runSingleSource(u);
        double sum2 = 0.0;
        G.forNodes([&](node v){
            if (u != v) {
                sum2 += ectd.runSinglePair(u,v);
            }
        });
        DEBUG("sum1 = ", sum1);
        DEBUG("sum2 = ", sum2);
        tlx::unused(sum1, sum2);
    //	INFO("Avg. relative error: ", error);
    }
}

} /* namespace NetworKit */
