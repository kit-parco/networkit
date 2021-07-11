// no-networkit-format
/*
 * MultiscaleBackboneGTest.cpp
 *
 *  Created on: 20.06.2014
 *      Author: Gerd Lindner
 */

#include <gtest/gtest.h>

#include <networkit/sparsification/Sparsifiers.hpp>
#include <networkit/sparsification/MultiscaleScore.hpp>


namespace NetworKit {

class MultiscaleBackboneGTest: public testing::Test {};

TEST_F(MultiscaleBackboneGTest, testSimpleMultiscaleBackbone) {
    Graph g(6, true, false);

    g.setWeight(0, 1, 1.0);
    g.setWeight(0, 2, 5.0);
    g.setWeight(0, 3, 10);
    g.setWeight(0, 4, 20.0);
    g.setWeight(4, 5, 1.0);
    g.setWeight(3, 5, 0.5);
    g.indexEdges();

    MultiscaleScore scorer(g, std::vector<double>());
    EXPECT_NEAR(0.9121756, scorer.getProbability(4, 0.5555), 1e-5) << "faulty probability calculation";
    /**
     * a01 = 0.91896
     * a02 = 0.639212
     * a03 = 0.376716
     * a04 = 0.0878244
     * a45 = 0.952381
     * a40 = 0.047619
     * a35 = 0.952381
     * a30 = 0.047619
     * a54 = 0.33333333
     * a53 = 0.66666666
     */

    //Compare the backbone graph to the expected backbone.

    MultiscaleSparsifier sparsifier(g, 0.5);
    sparsifier.run();
    Graph b = sparsifier.getGraph();

    EXPECT_EQ(3, b.numberOfEdges());
    EXPECT_TRUE(b.hasEdge(0, 4));
    EXPECT_TRUE(b.hasEdge(0, 3));
    EXPECT_TRUE(b.hasEdge(4, 5));

    MultiscaleSparsifier sparsifier2(g, 0.7);
    sparsifier2.run();
    b = sparsifier2.getGraph();
    EXPECT_EQ(2, b.numberOfEdges());
    EXPECT_TRUE(b.hasEdge(0, 4));
    EXPECT_TRUE(b.hasEdge(0, 3));
}

}
/* namespace NetworKit */
