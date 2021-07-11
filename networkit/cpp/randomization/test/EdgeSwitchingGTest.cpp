/*
 *  EdgeSwitchingGTest.cpp
 *
 *  Created on: Oct 20, 2019
 *  Author: Manuel Penschuck <networkit@manuel.jetzt>
 */

#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/randomization/EdgeSwitching.hpp>

namespace NetworKit {

class EdgeSwitchingGTest : public testing::TestWithParam<bool> {
protected:
    bool isDirected() const { return GetParam(); }

    /**
     * Creates a graph with two disjoint star graphs. In the directed case edges go from the
     * hubs to the satellites.
     *
     * Cluster 0: contains nodes [0, numberOfSats1], where node 0 ist the hub
     * Cluster 1: contains nodes [numberOfSats1+1, numberOfSats1+numberOfSats2+2], where node
     * numberOfSats1+1 is the hub
     */
    Graph createDualStarGraph(count numberOfSats1, count numberOfSats2) {
        Graph G(numberOfSats1 + numberOfSats2 + 2, false, isDirected());

        const node hub1 = 0;
        const node hub2 = numberOfSats1 + 1;

        for (node u = 1; u <= numberOfSats1; ++u)
            G.addEdge(hub1, hub1 + u);

        for (node u = 1; u <= numberOfSats2; ++u)
            G.addEdge(hub2, hub2 + u);

        return G;
    }

    /**
     * Create a circle graph with edges (i, i+1) and (n-1, 0).
     */
    Graph createCircle(count numberOfNodes) {
        Graph G(numberOfNodes, false, isDirected());

        G.addEdge(numberOfNodes - 1, 0);
        for (node u = 0; u < numberOfNodes - 1; ++u)
            G.addEdge(u, u + 1);

        return G;
    }
};

INSTANTIATE_TEST_SUITE_P(EdgeSwitchingDirected, EdgeSwitchingGTest, testing::Bool());

TEST_P(EdgeSwitchingGTest, testSatellites) {
    Aux::Random::setSeed(1, true);

    // Here we are switching two disjoint star graphs of equal size. As a result, expected half
    // of all switches are illegal. So we're checking proper rejection.
    constexpr count kNumSatellites = 100;
    constexpr count kNumSwitches = 100 * kNumSatellites;
    const node hub1 = 0, hub2 = kNumSatellites + 1;
    auto G = createDualStarGraph(kNumSatellites, kNumSatellites);

    // Check that hubs have correct degree and are unconnected
    ASSERT_FALSE(G.hasEdge(hub1, hub2));
    ASSERT_EQ(G.degree(hub1), kNumSatellites);
    ASSERT_EQ(G.degree(hub2), kNumSatellites);

    EdgeSwitchingInPlace algo(G, kNumSatellites); // carry out perturbation in-place
    ASSERT_EQ(algo.getNumberOfSwitchesPerEdge(), kNumSatellites);
    algo.setNumberOfSwitchesPerEdge(kNumSatellites / 2);
    ASSERT_EQ(algo.getNumberOfSwitchesPerEdge(), kNumSatellites / 2);

    for (int iter = 0; iter < 10; ++iter) {
        const auto numSwapsBefore = algo.getNumberOfAffectedEdges();
        algo.run();
        const auto numSuccessfulSwaps = (algo.getNumberOfAffectedEdges() - numSwapsBefore) / 2;

        // A directed switch is illegal with probability ~0.5 (if two edges are selected from the
        // same cluster). In case of an undirected swap, only one of the two possible resulting
        // topologies is illegal. So ~0.25 will be rejected.
        ASSERT_NEAR(numSuccessfulSwaps, kNumSwitches / (isDirected() ? 2 : 4),
                    kNumSwitches / (isDirected() ? 4 : 8));

        // Ensure a sufficient number of nodes changed their hub
        ASSERT_EQ(G.degree(0), kNumSatellites);
        auto range = G.neighborRange(0);
        const auto numSwitchedNodes =
            std::count_if(range.begin(), range.end(), [&](node u) { return u > kNumSatellites; });
        ASSERT_NEAR(numSwitchedNodes, kNumSatellites / 2, kNumSatellites / 5);
    }
}

TEST_P(EdgeSwitchingGTest, testCircle) {
    Aux::Random::setSeed(1, true);

    constexpr count kNumNodes = 100;
    constexpr count kNumSwitches = 100 * kNumNodes;
    auto Ginput = createCircle(kNumNodes);

    EdgeSwitching algo(Ginput, kNumNodes + 1, false);
    ASSERT_EQ(algo.getNumberOfSwitchesPerEdge(), kNumNodes + 1);
    algo.setNumberOfSwitchesPerEdge(kNumNodes);
    ASSERT_EQ(algo.getNumberOfSwitchesPerEdge(), kNumNodes);

    for (int iter = 0; iter < 10; ++iter) {
        const auto numSwapsBefore = algo.getNumberOfAffectedEdges();
        algo.run();
        const auto numSuccessfulSwaps = (algo.getNumberOfAffectedEdges() - numSwapsBefore) / 2;

        // Check that input graph is untouched
        Ginput.forNodes([&](node u) {
            ASSERT_TRUE(Ginput.hasEdge(u, (u + 1) % kNumNodes));
            ASSERT_TRUE(Ginput.hasEdge((u + kNumNodes - 1) % kNumNodes, u));
        });

        const auto &resultGraph = algo.getGraph();

        // For any edge there are exactly 2 illegal partners, so we expect an success rate of
        const auto expectedSuccessRate = static_cast<double>(kNumNodes - 2) / kNumNodes;
        const auto expectedSuccessfulSwaps = static_cast<count>(kNumSwitches * expectedSuccessRate);
        ASSERT_NEAR(numSuccessfulSwaps, expectedSuccessfulSwaps, expectedSuccessfulSwaps * 5 / 100);

        // In a truely random graph, the edge (u, u+1) exists with probability at most
        // 2.0 / kNumNodes, so we expect at most 2 such edges.
        count numContigousEdges = 0;
        resultGraph.forEdges([&](node u, node v) {
            if (u > v)
                std::swap(u, v);
            if (!u) {
                numContigousEdges += (v == kNumNodes - 1);
            } else {
                numContigousEdges += (u + 1 == v);
            }
        });

        resultGraph.forNodes([&](node u) {
            if (isDirected()) {
                ASSERT_EQ(resultGraph.degreeOut(u), 1);
                ASSERT_EQ(resultGraph.degreeIn(u), 1);
            } else {
                ASSERT_EQ(resultGraph.degree(u), 2);
            }
        });

        // Error prob < 1e-5
        ASSERT_LE(numContigousEdges, 10);
    }
}

TEST_P(EdgeSwitchingGTest, testPreprocessing) {
    Aux::Random::setSeed(1, true);

    constexpr count kNumNodes = 100;
    auto Ginput = createCircle(kNumNodes);

    for (int iter = 0; iter < 10; ++iter) {
        // We only use preprocessing as randomization. Since we have a two-regular graph,
        // preprocessing already yields a uniform sample.
        EdgeSwitching algo(Ginput, 0, true);
        ASSERT_EQ(algo.getNumberOfSwitchesPerEdge(), 0);

        algo.run();

        // In a truely random graph, the edge (u, u+1) exists with probability at most
        // 2.0 / kNumNodes, so we expect at most 2 such edges.
        count numContigousEdges = 0;
        algo.getGraph().forEdges([&](node u, node v) {
            if (u > v)
                std::swap(u, v);
            if (!u) {
                numContigousEdges += (v == kNumNodes - 1);
            } else {
                numContigousEdges += (u + 1 == v);
            }
        });

        ASSERT_LE(numContigousEdges, 10);
    }
}

} // namespace NetworKit
