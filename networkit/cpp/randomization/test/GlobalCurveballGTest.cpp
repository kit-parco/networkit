/*
 * GlobalCurveballGTest.cpp
 *
 *  Created on: 24.05.2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#include <gtest/gtest.h>

#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/randomization/GlobalCurveball.hpp>

namespace NetworKit {

class GlobalCurveballGTest : public ::testing::Test {
protected:
    void checkWithUndirectedGraph(Graph &);
    void checkWithDirectedGraph(Graph &G, bool selfLoops);
};

void GlobalCurveballGTest::checkWithUndirectedGraph(Graph &G) {
    ASSERT_FALSE(G.isDirected());

    node numNodes = G.numberOfNodes();
    const count numTrades = 5;

    std::vector<node> degrees(numNodes + 1);

    // Add edge to node 0, if isolated node
    // If 0 itself is isolated, add new node and connect 0 to it
    G.forNodes([&](node u) {
        if (G.degree(u) > 0)
            degrees[u] = G.degree(u);
        else {
            if (u == 0) {
                numNodes++;
                G.addEdge(0, numNodes - 1);
                degrees[0]++;
                degrees[numNodes - 1] = 1;
            } else {
                G.addEdge(u, 0);
                degrees[0]++;
                degrees[u] = 1;
            }
        }
    });

    GlobalCurveball algo(G, numTrades);
    algo.run();

    // check degrees
    Graph outG = algo.getGraph();
    outG.forNodes([&](node u) { ASSERT_EQ(degrees[u], outG.degree(u)); });
}

void GlobalCurveballGTest::checkWithDirectedGraph(Graph &G, bool selfLoops) {
    ASSERT_TRUE(G.isDirected());

    node numNodes = G.numberOfNodes();
    const count numTrades = 5;

    std::vector<node> degreesIn(numNodes + 1);
    std::vector<node> degreesOut(numNodes + 1);

    // Add edge to node 0, if isolated node
    // If 0 itself is isolated, add new node and connect 0 to it
    G.forNodes([&](node u) {
        if (G.degreeIn(u) > 0 || G.degreeOut(u) > 0) {
            degreesIn[u] = G.degreeIn(u);
            degreesOut[u] = G.degreeOut(u);
        } else {
            if (u == 0) {
                numNodes++;
                G.addEdge(0, numNodes - 1);
                degreesOut[0]++;
                degreesIn[numNodes - 1] = 1;
            } else {
                G.addEdge(u, 0);
                degreesIn[0]++;
                degreesOut[u] = 1;
            }
        }
    });

    GlobalCurveball algo(G, numTrades, selfLoops);
    algo.run();

    // check degrees
    Graph outG = algo.getGraph();
    outG.forNodes([&](node u) {
        ASSERT_EQ(degreesIn[u], outG.degreeIn(u));
        ASSERT_EQ(degreesOut[u], outG.degreeOut(u));
    });

    if (selfLoops) {
        ASSERT_GT(outG.numberOfSelfLoops(), 0);
    } else {
        ASSERT_EQ(outG.numberOfSelfLoops(), 0);
    }
}

TEST_F(GlobalCurveballGTest, testCurveballUndirectedErdosRenyi) {
    Aux::Random::setSeed(1, false);

    node numNodes = 1000;
    ErdosRenyiGenerator generator(numNodes, 0.01);
    Graph G = generator.generate();

    this->checkWithUndirectedGraph(G);
}

TEST_F(GlobalCurveballGTest, testCurveballDirectedErdosRenyi) {
    Aux::Random::setSeed(1, false);

    node numNodes = 1000;
    ErdosRenyiGenerator generator(numNodes, 0.01, true, false);
    Graph G = generator.generate();

    this->checkWithDirectedGraph(G, false);
}

TEST_F(GlobalCurveballGTest, testCurveballDirectedErdosRenyiSelfLoops) {
    Aux::Random::setSeed(1, false);

    node numNodes = 1000;
    ErdosRenyiGenerator generator(numNodes, 0.01, true, true);
    Graph G = generator.generate();

    this->checkWithDirectedGraph(G, true);
}

TEST_F(GlobalCurveballGTest, testCurveballHyperbolic) {
    Aux::Random::setSeed(1, false);

    node numNodes = 1000;
    HyperbolicGenerator generator(numNodes);
    Graph G = generator.generate();

    this->checkWithUndirectedGraph(G);
}

} // namespace NetworKit
