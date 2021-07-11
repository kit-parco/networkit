// no-networkit-format
/*
 * AlgebraicBFSGTest.cpp
 *
 *  Created on: Jun 7, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <gtest/gtest.h>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/algebraic/algorithms/AlgebraicBFS.hpp>
#include <networkit/io/METISGraphReader.hpp>

#include <networkit/auxiliary/Timer.hpp>

namespace NetworKit {

class AlgebraicBFSGTest : public testing::Test {};

TEST(AlgebraicBFSGTest, testOnToyGraph) {
    Graph G(5, false, true);
    G.addEdge(0, 1);
    G.addEdge(0, 3);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(1, 4);
    G.addEdge(2, 1);
    G.addEdge(3, 2);
    G.addEdge(3, 4);
    G.addEdge(4, 0);
    G.addEdge(4, 2);

    AlgebraicBFS<CSRMatrix> bfs(G, 0);
    bfs.run();

    EXPECT_EQ(0, bfs.distance(0));
    EXPECT_EQ(1, bfs.distance(1));
    EXPECT_EQ(2, bfs.distance(2));
    EXPECT_EQ(1, bfs.distance(3));
    EXPECT_EQ(2, bfs.distance(4));

    G = Graph(7, false, true);
    G.addEdge(0,1);
    G.addEdge(0,3);
    G.addEdge(1,4);
    G.addEdge(1,6);
    G.addEdge(2,5);
    G.addEdge(3,0);
    G.addEdge(3,2);
    G.addEdge(4,5);
    G.addEdge(5,2);
    G.addEdge(6,2);
    G.addEdge(6,3);

    bfs = AlgebraicBFS<CSRMatrix>(G,3);
    bfs.run();

    EXPECT_EQ(1, bfs.distance(0));
    EXPECT_EQ(2, bfs.distance(1));
    EXPECT_EQ(1, bfs.distance(2));
    EXPECT_EQ(0, bfs.distance(3));
    EXPECT_EQ(3, bfs.distance(4));
    EXPECT_EQ(2, bfs.distance(5));
    EXPECT_EQ(3, bfs.distance(6));
}

TEST(AlgebraicBFSGTest, benchmarkBFS) {
    METISGraphReader reader;
    Graph G = reader.read("input/caidaRouterLevel.graph");

    Aux::Timer t;
    BFS bfs(G, 0, false);
    t.start();
    bfs.run();
    t.stop();
    INFO("Bfs took ", t.elapsedMilliseconds());

    AlgebraicBFS<CSRMatrix> algebraicBfs(G, 0);

    t.start();
    algebraicBfs.run();
    t.stop();
    INFO("AlgebraicBFS took ", t.elapsedMilliseconds());

    G.forNodes([&](node u) {
        if (algebraicBfs.distance(u) == std::numeric_limits<double>::infinity()) {
            EXPECT_EQ(bfs.distance(u), std::numeric_limits<double>::max());
        } else {
            EXPECT_EQ(bfs.distance(u), algebraicBfs.distance(u));
        }
    });

}

} /* namespace NetworKit */
