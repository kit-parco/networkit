// no-networkit-format
/*
 * BellmanFordGTest.cpp
 *
 *  Created on: Jun 6, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <gtest/gtest.h>

#include <networkit/io/METISGraphReader.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/graph/Graph.hpp>

#include <networkit/algebraic/CSRMatrix.hpp>

#include <networkit/distance/Dijkstra.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>
#include <networkit/algebraic/algorithms/AlgebraicBellmanFord.hpp>

namespace NetworKit {

class AlgebraicBellmanFordGTest : public testing::Test {
public:
    AlgebraicBellmanFordGTest() = default;
    virtual ~AlgebraicBellmanFordGTest() = default;

protected:
    std::vector<double> classicBF(const Graph& graph, node s) const;
};

TEST_F(AlgebraicBellmanFordGTest, testOnToyGraph) {
    Graph G(5, true, true);
    G.addEdge(0, 1, 6);
    G.addEdge(0, 3, 7);
    G.addEdge(1, 2, 5);
    G.addEdge(1, 3, 8);
    G.addEdge(1, 4, -4);
    G.addEdge(2, 1, -2);
    G.addEdge(3, 2, -3);
    G.addEdge(3, 4, 9);
    G.addEdge(4, 0, 2);
    G.addEdge(4, 2, 7);

    AlgebraicBellmanFord<CSRMatrix> bf(G, 0);

    bf.run();

    EXPECT_EQ(0, bf.distance(0));
    EXPECT_EQ(2, bf.distance(1));
    EXPECT_EQ(4, bf.distance(2));
    EXPECT_EQ(7, bf.distance(3));
    EXPECT_EQ(-2, bf.distance(4));

    EXPECT_FALSE(bf.hasNegativeCycle());
}

TEST_F(AlgebraicBellmanFordGTest, benchmark) {
    METISGraphReader reader;
    Graph graph = reader.read("input/PGPgiantcompo.graph");

    AlgebraicBellmanFord<CSRMatrix> bf(graph, 0);

    Aux::Timer t;

    t.start();
    std::vector<double> classicResult = classicBF(graph, 0);
    t.stop();

    INFO("Classic BF took ", t.elapsedMilliseconds(), " ms");

    t.start();
    bf.run();
    t.stop();

    INFO("Algebraic BF took ", t.elapsedMilliseconds(), " ms.");

    Dijkstra dijkstra(graph, 0);

    t.start();
    dijkstra.run();
    t.stop();
    INFO("Dijkstra took ", t.elapsedMilliseconds(), " ms.");

    for (index i = 0; i < classicResult.size(); ++i) {
        EXPECT_EQ(dijkstra.distance(i), bf.distance(i));
        EXPECT_EQ(classicResult[i], bf.distance(i));
    }
}

std::vector<double> AlgebraicBellmanFordGTest::classicBF(const Graph& graph, node s) const {
    std::vector<double> dist(graph.numberOfNodes(), std::numeric_limits<double>::infinity());
    dist[s] = 0;

    for (index i = 1; i < graph.numberOfNodes(); ++i) {
        graph.forNodes([&](node u) {
            if (dist[u] != std::numeric_limits<double>::infinity()) {
                graph.forNeighborsOf(u, [&](node v, edgeweight w) {
                    dist[v] = std::min(dist[v], dist[u]+w);
                });
            }
        });
    }

    return dist;
}

} /* namespace NetworKit */
