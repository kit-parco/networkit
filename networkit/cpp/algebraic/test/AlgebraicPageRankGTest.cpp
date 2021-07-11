// no-networkit-format
/*
 * AlgebraicPageRankGTest.cpp
 *
 *  Created on: Jun 20, 2016
 *      Author: Michael
 */

#include <gtest/gtest.h>

#include <networkit/algebraic/algorithms/AlgebraicPageRank.hpp>
#include <networkit/io/SNAPGraphReader.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/auxiliary/Timer.hpp>

#include <networkit/centrality/PageRank.hpp>

namespace NetworKit {

class AlgebraicPageRankGTest : public testing::Test {};

TEST(AlgebraicPageRankGTest, testPageRankDirected) {
    SNAPGraphReader reader;
    Graph G = reader.read("input/wiki-Vote.txt");

    Aux::Timer t;
    t.start();
    AlgebraicPageRank<CSRMatrix> apr(G);
    t.stop();

    INFO("Initializing AlgebraicPageRank took ", t.elapsedMicroseconds() / 1000.0);

    t.start();
    apr.run();
    t.stop();
    INFO("Computing page rank algebraically took ", t.elapsedMicroseconds() / 1000.0);
    std::vector<std::pair<node, double> > pr_ranking = apr.ranking();

    const double tol = 1e-3;
    EXPECT_EQ(pr_ranking[0].first, 699u);
    EXPECT_NEAR(pr_ranking[0].second, 0.00432, tol);

    /*t.start();
    PageRank pr(G);
    t.stop();
    INFO("Initilaizing graph-theoretic PageRank took ", t.elapsedMicroseconds() / 1000.0);

    t.start();
    pr.run();
    t.stop();

    INFO("Computing graph-theoretic page rank took ", t.elapsedMicroseconds() / 1000.0);*/
}

TEST(AlgebraicPageRankGTest, testPageRankCentrality) {
 /* Graph:
    0    3   6
     \  / \ /
      2 -- 5
     /  \ / \
    1    4   7

    Edges in the upper row have weight 3,
    the edge in the middle row has weight 1.5,
    edges in the lower row have weight 2.
 */
    count n = 8;
    Graph G(n, true);

    G.addEdge(0, 2, 3);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 3, 3);
    G.addEdge(2, 4, 2);
    G.addEdge(2, 5, 1.5);
    G.addEdge(3, 5, 3);
    G.addEdge(4, 5, 2);
    G.addEdge(5, 6, 3);
    G.addEdge(5, 7, 2);

    double damp = 0.85;
    AlgebraicPageRank<CSRMatrix> apr(G, damp);
    apr.run();
    std::vector<double> cen = apr.scores();

    // compare to Matlab results
    const double tol = 1e-4;
    EXPECT_NEAR(0.0753, fabs(cen[0]), tol);
    EXPECT_NEAR(0.0565, fabs(cen[1]), tol);
    EXPECT_NEAR(0.2552, fabs(cen[2]), tol);
    EXPECT_NEAR(0.1319, fabs(cen[3]), tol);
    EXPECT_NEAR(0.0942, fabs(cen[4]), tol);
    EXPECT_NEAR(0.2552, fabs(cen[5]), tol);
    EXPECT_NEAR(0.0753, fabs(cen[6]), tol);
    EXPECT_NEAR(0.0565, fabs(cen[7]), tol);
}

} /* namespace NetworKit */
