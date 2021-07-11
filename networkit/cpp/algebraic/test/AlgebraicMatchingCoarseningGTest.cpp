// no-networkit-format
/*
 * AlgebraicMatchingCoarseningGTest.cpp
 *
 *  Created on: Jul 12, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/algebraic/algorithms/AlgebraicMatchingCoarsening.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/matching/LocalMaxMatcher.hpp>
#include <networkit/coarsening/MatchingCoarsening.hpp>

namespace NetworKit {

class AlgebraicMatchingCoarseningGTest : public testing::Test{};


TEST_F(AlgebraicMatchingCoarseningGTest, testContraction) {
    METISGraphReader reader;
    Graph G = reader.read("input/wing.graph");

    LocalMaxMatcher matcher(G);
    matcher.run();
    Matching matching = matcher.getMatching();
    ASSERT_TRUE(matching.isProper(G));

    Aux::Timer t;
    t.start();
    AlgebraicMatchingCoarsening<CSRMatrix> amc(G, matching);
    amc.run();
    t.stop();
    INFO("Algebraic matching coarsening took ", t.elapsedTag());
    Graph coarseG = amc.getCoarseGraph();
    std::vector<node> amcFineToCoarse = amc.getFineToCoarseNodeMapping();

    t.start();
    MatchingCoarsening mc(G, matching);
    mc.run();
    t.stop();
    INFO("Graph theoretic matching coarsening took ", t.elapsedTag());
    std::vector<node> mcFineToCoarse = mc.getFineToCoarseNodeMapping();

    G.forNodes([&](node u) {
        EXPECT_EQ(mcFineToCoarse[u], amcFineToCoarse[u]);
    });

    EXPECT_GE(coarseG.numberOfNodes(), 0.5*G.numberOfNodes());
    EXPECT_EQ(G.totalEdgeWeight(), coarseG.totalEdgeWeight());
    if (G.numberOfEdges() > 0) {
        EXPECT_LT(coarseG.numberOfNodes(), G.numberOfNodes());
    }
}

} /* namespace NetworKit */
