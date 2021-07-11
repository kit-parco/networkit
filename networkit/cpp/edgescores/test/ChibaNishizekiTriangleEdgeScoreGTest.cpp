// no-networkit-format
/*
 * ChibaNishizekiEdgeScoreGTest.cpp
 *
 *  Created on: 23.05.2014
 *      Author: Gerd Lindner
 */

#include <gtest/gtest.h>

#include <networkit/edgescores/ChibaNishizekiTriangleEdgeScore.hpp>
#include <networkit/edgescores/TriangleEdgeScore.hpp>

namespace NetworKit {

class ChibaNishizekiTriangleEdgeScoreGTest: public testing::Test {};

TEST_F(ChibaNishizekiTriangleEdgeScoreGTest, testTriangleCountsTrivial) {
    Graph g(5);

    g.addEdge(0,1);
    g.addEdge(0,2);
    g.addEdge(1,2);

    g.indexEdges();

    ChibaNishizekiTriangleEdgeScore counter(g);
    counter.run();
    std::vector<count> counts = counter.scores();

    EXPECT_EQ(1, (counts[g.edgeId(0,1)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(0,2)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(1,2)])) << "wrong triangle count";
}

TEST_F(ChibaNishizekiTriangleEdgeScoreGTest, testNewTriangleCountsTrivial) {
    Graph g(5);

    g.addEdge(0,1);
    g.addEdge(0,2);
    g.addEdge(1,2);

    g.indexEdges();

    TriangleEdgeScore counter(g);
    counter.run();
    std::vector<count> counts = counter.scores();

    EXPECT_EQ(1, (counts[g.edgeId(0,1)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(0,2)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(1,2)])) << "wrong triangle count";
    //EXPECT_EQ(0, (counts[g.edgeId(2,3)])) << "wrong triangle count";
    //TODO: edge ids for non-existing edges currently result in unexpected behaviour.
}


TEST_F(ChibaNishizekiTriangleEdgeScoreGTest, testTriangleCountsSimple) {
    int64_t n = 6;
    Graph g(n);

    g.addEdge(0,1);
    g.addEdge(0,2);
    g.addEdge(1,2);

    g.addEdge(0,4);
    g.addEdge(0,3);
    g.addEdge(3,4);

    g.addEdge(0,5);
    g.addEdge(4,5);

    g.indexEdges();

    EXPECT_EQ(8, g.numberOfEdges()) << "wrong edge count";

    ChibaNishizekiTriangleEdgeScore counter(g);
    counter.run();
    std::vector<count> counts = counter.scores();

    EXPECT_EQ(6, g.numberOfNodes()) << "undesired side effect";
    EXPECT_EQ(8, g.numberOfEdges()) << "undesired side effect";

    EXPECT_EQ(8, counts.size()) << "wrong triangle count map size";
    EXPECT_EQ(1, (counts[g.edgeId(0,1)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(0,2)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(1,2)])) << "wrong triangle count";

    EXPECT_EQ(1, (counts[g.edgeId(0,3)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(3,4)])) << "wrong triangle count";

    EXPECT_EQ(2, (counts[g.edgeId(0,4)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(0,5)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(4,5)])) << "wrong triangle count";

    EXPECT_EQ(1, (counts[g.edgeId(1,0)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(2,0)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(2,1)])) << "wrong triangle count";

    EXPECT_EQ(1, (counts[g.edgeId(3,0)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(4,3)])) << "wrong triangle count";

    EXPECT_EQ(2, (counts[g.edgeId(4,0)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(5,0)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(5,4)])) << "wrong triangle count";
}

TEST_F(ChibaNishizekiTriangleEdgeScoreGTest, testNewTriangleCountsSimple) {
    int64_t n = 6;
    Graph g(n);

    g.addEdge(0,1);
    g.addEdge(0,2);
    g.addEdge(1,2);

    g.addEdge(0,4);
    g.addEdge(0,3);
    g.addEdge(3,4);

    g.addEdge(0,5);
    g.addEdge(4,5);

    g.indexEdges();

    EXPECT_EQ(8, g.numberOfEdges()) << "wrong edge count";

    TriangleEdgeScore counter(g);
    counter.run();
    std::vector<count> counts = counter.scores();

    EXPECT_EQ(6, g.numberOfNodes()) << "undesired side effect";
    EXPECT_EQ(8, g.numberOfEdges()) << "undesired side effect";

    EXPECT_EQ(8, counts.size()) << "wrong triangle count map size";
    EXPECT_EQ(1, (counts[g.edgeId(0,1)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(0,2)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(1,2)])) << "wrong triangle count";

    EXPECT_EQ(1, (counts[g.edgeId(0,3)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(3,4)])) << "wrong triangle count";

    EXPECT_EQ(2, (counts[g.edgeId(0,4)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(0,5)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(4,5)])) << "wrong triangle count";

    EXPECT_EQ(1, (counts[g.edgeId(1,0)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(2,0)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(2,1)])) << "wrong triangle count";

    EXPECT_EQ(1, (counts[g.edgeId(3,0)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(4,3)])) << "wrong triangle count";

    EXPECT_EQ(2, (counts[g.edgeId(4,0)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(5,0)])) << "wrong triangle count";
    EXPECT_EQ(1, (counts[g.edgeId(5,4)])) << "wrong triangle count";
}

}

/* namespace NetworKit */
