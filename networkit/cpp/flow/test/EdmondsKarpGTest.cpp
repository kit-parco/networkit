// no-networkit-format
/*
* EdmondsKarpGTest.cpp
 *
 *  Created on: 13.06.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <gtest/gtest.h>

#include <networkit/flow/EdmondsKarp.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class EdmondsKarpGTest : public testing::Test{};

TEST_F(EdmondsKarpGTest, testEdmondsKarpP1) {
    Graph G(7, false);
    G.addEdge(0,1);
    G.addEdge(0,2);
    G.addEdge(0,3);
    G.addEdge(1,2);
    G.addEdge(1,4);
    G.addEdge(2,3);
    G.addEdge(2,4);
    G.addEdge(3,4);
    G.addEdge(3,5);
    G.addEdge(4,6);
    G.addEdge(5,6);

    G.indexEdges();

    EdmondsKarp edKa(G, 0, 6);
    edKa.run();
    EXPECT_EQ(2, edKa.getMaxFlow()) << "max flow is not correct";

    EXPECT_EQ(1, edKa.getFlow(4, 6));
    EXPECT_EQ(1, edKa.getFlow(5, 6));

    std::vector<node> sourceSet(edKa.getSourceSet());

    EXPECT_NE(std::find(sourceSet.begin(), sourceSet.end(), 0), sourceSet.end());
    EXPECT_NE(std::find(sourceSet.begin(), sourceSet.end(), 1), sourceSet.end());
    EXPECT_NE(std::find(sourceSet.begin(), sourceSet.end(), 2), sourceSet.end());
    EXPECT_NE(std::find(sourceSet.begin(), sourceSet.end(), 3), sourceSet.end());
    EXPECT_NE(std::find(sourceSet.begin(), sourceSet.end(), 4), sourceSet.end());

    EXPECT_EQ(std::find(sourceSet.begin(), sourceSet.end(), 5), sourceSet.end());
    EXPECT_EQ(std::find(sourceSet.begin(), sourceSet.end(), 6), sourceSet.end());
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpTwoPaths) {
    Graph G(11);

    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.addEdge(1, 4);
    G.addEdge(1, 5);
    G.addEdge(5, 6);
    G.addEdge(6, 7);
    G.addEdge(7, 8);
    G.addEdge(8, 9);
    G.addEdge(4, 10);
    G.addEdge(9, 10);

    G.indexEdges();

    EdmondsKarp edKa(G, 0, 10);
    edKa.run();

    EXPECT_EQ(2, edKa.getMaxFlow());
    EXPECT_EQ(0, edKa.getFlow(1, 4));
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpP2) {
    Graph G(6, true);
    G.addEdge(0,1, 5);
    G.addEdge(0,2, 15);
    G.addEdge(1,3, 5);
    G.addEdge(1,4, 5);
    G.addEdge(2,3, 5);
    G.addEdge(2, 4, 5);
    G.addEdge(3,5, 15);
    G.addEdge(4,5, 5);

    G.indexEdges();

    EdmondsKarp edKa(G, 0, 5);
    edKa.run();

    EXPECT_EQ(15, edKa.getMaxFlow()) << "max flow is not correct";
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpUnconnected) {
    Graph G(6, true);
    G.addEdge(0,1, 5);
    G.addEdge(0,2, 15);
    G.addEdge(1,2, 5);
    G.addEdge(3, 4, 5);
    G.addEdge(3,5, 15);
    G.addEdge(4,5, 5);

    G.indexEdges();

    EdmondsKarp edKa(G, 0, 5);
    edKa.run();
    EXPECT_EQ(0, edKa.getMaxFlow()) << "max flow is not correct";
}

} /* namespace NetworKit */
