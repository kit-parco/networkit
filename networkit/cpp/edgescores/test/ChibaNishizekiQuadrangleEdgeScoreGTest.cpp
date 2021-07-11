// no-networkit-format
/*
 * ChibaNishizekiEdgeScoreGTest.cpp
 *
 *  Created on: 23.05.2014
 *      Author: Gerd Lindner
 */

#include <gtest/gtest.h>

#include <networkit/edgescores/ChibaNishizekiQuadrangleEdgeScore.hpp>

namespace NetworKit {

class ChibaNishizekiQuadrangleEdgeScoreGTest: public testing::Test {};

TEST_F(ChibaNishizekiQuadrangleEdgeScoreGTest, testQuadrangleCountsTrivial) {
    Graph g(5);

    g.addEdge(0,1);
    g.addEdge(0,2);
    g.addEdge(1,3);
    g.addEdge(2,3);

    g.indexEdges();

    ChibaNishizekiQuadrangleEdgeScore counter(g);
    counter.run();
    std::vector<count> counts = counter.scores();

    EXPECT_EQ(1, (counts[g.edgeId(0,1)])) << "wrong quadrangle count";
    EXPECT_EQ(1, (counts[g.edgeId(0,2)])) << "wrong quadrangle count";
    EXPECT_EQ(1, (counts[g.edgeId(1,3)])) << "wrong quadrangle count";
    EXPECT_EQ(1, (counts[g.edgeId(2,3)])) << "wrong quadrangle count";
    //EXPECT_EQ(0, (counts[g.edgeId(2,3)])) << "wrong quadrangle count";
    //TODO: edge ids for non-existing edges currently result in unexpected behaviour.
}

TEST_F(ChibaNishizekiQuadrangleEdgeScoreGTest, testQuadrangleCountsSimple) {
    count n = 7;
    Graph g(n);

    g.addEdge(0,1);
    g.addEdge(0,3);
    g.addEdge(0,4);
    g.addEdge(0,6);
    g.addEdge(1,2);
    g.addEdge(1,3);
    g.addEdge(2,3);
    g.addEdge(3,5);
    g.addEdge(3,6);
    g.addEdge(4,5);

    g.indexEdges();

    EXPECT_EQ(10, g.numberOfEdges()) << "wrong edge count";

    ChibaNishizekiQuadrangleEdgeScore counter(g);
    counter.run();
    std::vector<count> counts = counter.scores();

    EXPECT_EQ(7, g.numberOfNodes()) << "undesired side effect";
    EXPECT_EQ(10, g.numberOfEdges()) << "undesired side effect";

    EXPECT_EQ(10, counts.size()) << "wrong quadrangle count map size";
    EXPECT_EQ(2, (counts[g.edgeId(0,1)])) << "wrong quadrangle count";
    EXPECT_EQ(2, (counts[g.edgeId(0,3)])) << "wrong quadrangle count";
    EXPECT_EQ(1, (counts[g.edgeId(0,4)])) << "wrong quadrangle count";
    EXPECT_EQ(1, (counts[g.edgeId(0,6)])) << "wrong quadrangle count";
    
    EXPECT_EQ(1, (counts[g.edgeId(1,2)])) << "wrong quadrangle count";
    EXPECT_EQ(1, (counts[g.edgeId(1,3)])) << "wrong quadrangle count";
    EXPECT_EQ(1, (counts[g.edgeId(2,3)])) << "wrong quadrangle count";

    EXPECT_EQ(1, (counts[g.edgeId(3,5)])) << "wrong quadrangle count";
    EXPECT_EQ(1, (counts[g.edgeId(3,6)])) << "wrong quadrangle count";
    EXPECT_EQ(1, (counts[g.edgeId(4,5)])) << "wrong quadrangle count";
}


}
/* namespace NetworKit */
