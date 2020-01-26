/*
 * LowToHighDirectedGraphGTest.cpp
 *
 * Created: 2020-01-12
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>

#include <algorithm>

#include <networkit/structures/LowToHighDirectedGraph.hpp>

namespace NetworKit {

class LowToHighDirectedGraphGTest : public testing::Test {
};

TEST_F(LowToHighDirectedGraphGTest, testForEdgesOf) {
    Graph G(5, true);
    G.addEdge(0, 1, 0.9);
    G.addEdge(0, 2, 0.8);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(2, 3);

    LowToHighDirectedGraph directedGraph(G);

    std::vector<node> neighbors;
    edgeweight weightedDegree = 0;
    directedGraph.forEdgesOf(0, [&](node, node v, edgeweight weight) {
        neighbors.push_back(v);
        weightedDegree += weight;
    });
    ASSERT_EQ(neighbors, std::vector<node>({1, 2}));
    ASSERT_NEAR(weightedDegree, 1.7, 1e-6);
}

TEST_F(LowToHighDirectedGraphGTest, testNeighborsOf) {
    Graph G(5);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 4);

    LowToHighDirectedGraph directedGraph(G);

    auto sortedNeighbors = [&](node u) {
        auto neighbors = directedGraph.neighborsOf(u);
        std::sort(neighbors.begin(), neighbors.end());
        return neighbors;
    };
    ASSERT_EQ(sortedNeighbors(0), std::vector<node>{2});
    ASSERT_EQ(sortedNeighbors(1), std::vector<node>{2});
    ASSERT_EQ(sortedNeighbors(2), std::vector<node>{});
    ASSERT_EQ(sortedNeighbors(3), std::vector<node>{2});
    ASSERT_EQ(sortedNeighbors(4), std::vector<node>{3});
}

}