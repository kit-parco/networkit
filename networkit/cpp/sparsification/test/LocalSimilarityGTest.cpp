// no-networkit-format
/*
 * LocalSimilarityGTest.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#include <gtest/gtest.h>

#include <networkit/sparsification/LocalSimilarityScore.hpp>
#include <networkit/edgescores/ChibaNishizekiTriangleEdgeScore.hpp>

namespace NetworKit {

class LocalSimilarityGTest: public testing::Test {};

TEST_F(LocalSimilarityGTest, testAttributeSimple) {
    Graph g(4);

    g.addEdge(0, 1);
    g.addEdge(0, 3);
    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.indexEdges();

    ChibaNishizekiTriangleEdgeScore triangleEdgeScore(g);
    triangleEdgeScore.run();
    std::vector<count> triangles = triangleEdgeScore.scores();

    LocalSimilarityScore localSim(g, triangles);
    localSim.run();
    std::vector<double> exp = localSim.scores();

    EXPECT_DOUBLE_EQ(1.0, exp[g.edgeId(0, 1)]);
    EXPECT_NEAR(0.36907025, exp[g.edgeId(0, 2)], 1e-7);
    EXPECT_DOUBLE_EQ(1.0, exp[g.edgeId(0, 3)]);
    EXPECT_DOUBLE_EQ(1.0, exp[g.edgeId(1, 2)]);
}

}
/* namespace NetworKit */
