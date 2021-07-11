// no-networkit-format
/*
 * SpanningEdgeCentralityGTest.cpp
 *
 *  Created on: Jan 17, 2016
 *      Author: Michael
 */

#include <gtest/gtest.h>

#include <networkit/centrality/SpanningEdgeCentrality.hpp>

#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class SpanningEdgeCentralityGTest : public testing::Test {};

TEST_F(SpanningEdgeCentralityGTest, testOnToyGraph) {
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
}

} /* namespace NetworKit */
