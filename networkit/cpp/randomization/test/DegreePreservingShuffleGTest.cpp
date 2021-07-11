/*
 * DegreePreservingShuffleGTest.cpp
 *
 *  Created on: 26.08.2018
 *      Author:  Manuel Penschuck <networkit@manuel.jetzt>
 */

#include <gtest/gtest.h>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/randomization/DegreePreservingShuffle.hpp>

namespace NetworKit {

class DegreePreservingShuffleGTest : public ::testing::Test {};

namespace CurveballDetails {

TEST_F(DegreePreservingShuffleGTest, UndirectedNoChange) {
    const index n = 10;

    Graph G(2 * n);

    // build clique between nodes [n, 2n)
    for (index i = n + 1; i < 2 * n; i++) {
        for (index j = n; j < i; j++) {
            G.addEdge(i, j);
        }
    }

    // Connect nodes i in [0; n) to i-many neighbours in the clique
    for (index i = 1; i < n; i++) {
        for (index d = 0; d < i; d++)
            G.addEdge(i, 2 * n - d - 2);
    }

    DegreePreservingShuffle dps(G);
    dps.run();

    const auto perm = dps.getPermutation();
    const auto Gperm = dps.getGraph();

    for (index i = 0; i < n - 1; i++) {
        ASSERT_EQ(G.degree(i), i);
        ASSERT_EQ(perm[i], i);
        ASSERT_EQ(Gperm.degree(i), i);
    }
}

TEST_F(DegreePreservingShuffleGTest, UndirectedErdosRenyi) {
    // we need those rather large graphs to test the parallel algo
    for (auto n : {50, 100, 500, 1000, 200000}) {
        auto G = ErdosRenyiGenerator(n, 5.0 / n).generate();

        DegreePreservingShuffle dps(G);
        dps.run();

        const auto perm = dps.getPermutation();
        const auto Gperm = dps.getGraph();

        G.forNodes([&](node u) {
            ASSERT_EQ(G.degree(u), G.degree(perm[u]));
            ASSERT_EQ(G.degree(u), Gperm.degree(u));
        });
    }
}

TEST_F(DegreePreservingShuffleGTest, DirectedErdosRenyi) {
    // we need those rather large graphs to test the parallel algo
    for (auto n : {50, 100, 500, 1000, 200000}) {
        auto G = ErdosRenyiGenerator(n, 5.0 / n, true, true).generate();

        DegreePreservingShuffle dps(G);
        dps.run();

        const auto perm = dps.getPermutation();
        const auto Gperm = dps.getGraph();

        G.forNodes([&](node u) {
            ASSERT_EQ(G.degreeIn(u), G.degreeIn(perm[u]));
            ASSERT_EQ(G.degreeIn(u), Gperm.degreeIn(u));
            ASSERT_EQ(G.degreeOut(u), G.degreeOut(perm[u]));
            ASSERT_EQ(G.degreeOut(u), Gperm.degreeOut(u));
        });
    }
}

} // namespace CurveballDetails
} // namespace NetworKit
