// no-networkit-format
/*
 * DynBetweennessGTest.cpp
 *
 *  Created on: 05.08.2014
 *      Author: ebergamini, cls
 */

#include <gtest/gtest.h>

#include <networkit/centrality/Betweenness.hpp>
#include <networkit/centrality/DynApproxBetweenness.hpp>
#include <networkit/centrality/ApproxBetweenness.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/centrality/DynBetweenness.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class DynBetweennessGTest: public testing::Test {};

TEST_F(DynBetweennessGTest, runDynApproxBetweennessSmallGraph) {
/* Graph:
0    3   6
    \  / \ /
    2    5
    /  \ / \
1    4   7
*/
    int n = 8;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(5, 7);

    //double epsilon = 0.01; // error
    double epsilon = 0.1; // error
    double delta = 0.1; // confidence
    DynApproxBetweenness dynbc(G, epsilon, delta);
    Betweenness bc(G);
    dynbc.run();
    bc.run();
    std::vector<double> dynbc_scores = dynbc.scores();
    std::vector<double> bc_scores = bc.scores();
    for(int i=0; i<n; i++) {
        DEBUG("Difference ", dynbc_scores[i]-bc_scores[i]/double(n*(n-1)));
    }
    std::vector<GraphEvent> batch;
    batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, 0, 6, 1.0));
    G.addEdge(batch[0].u, batch[0].v);
    bc.run();
    dynbc.updateBatch(batch);
    dynbc_scores = dynbc.scores();
    bc_scores = bc.scores();
    for(int i=0; i<n; i++) {
        DEBUG("Difference ", dynbc_scores[i]-bc_scores[i]/double(n*(n-1)));
    }

}


TEST_F(DynBetweennessGTest, runDynVsStatic) {
    METISGraphReader reader;
    //Graph G = reader.read("input/PGPgiantcompo.graph");
    Graph G = reader.read("input/celegans_metabolic.graph");
    count n = G.upperNodeIdBound();

    double epsilon = 0.1; // error
    double delta = 0.1; // confidence
    DEBUG("Initializing DynApproxBetweenness");
    DynApproxBetweenness dynbc(G, epsilon, delta, false);
    DEBUG("Initializing ApproxBetweenness");
    ApproxBetweenness bc(G, epsilon, delta);
    DEBUG("Running DynApproxBetweenness");
    dynbc.run();
    DEBUG("Running ApproxBetweenness");
    bc.run();
    std::vector<double> dynbc_scores = dynbc.scores();
    std::vector<double> bc_scores = bc.scores();
    double err1=0;
    for(count i=0; i<n; i++) {
        double x = dynbc_scores[i]-bc_scores[i];
        if (x > err1)
            err1 = x;
    }
    DEBUG("Before the edge insertion: ");
    std::vector<GraphEvent> batch;
    count nInsertions = 10, i = 0;
    while (i < nInsertions) {
        node v1 = GraphTools::randomNode(G);
        node v2 = GraphTools::randomNode(G);
        if (v1 != v2 && !G.hasEdge(v1, v2)) {
            G.addEdge(v1, v2);
            batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0));
            i++;
        }
    }
    DEBUG("Running ApproxBetweenness (again)");
    bc.run();
    DEBUG("Updating DynApproxBetweenness");
    dynbc.updateBatch(batch);
    DEBUG("Calling DynApproxBetweenness Scores");
    dynbc_scores = dynbc.scores();
    DEBUG("Calling ApproxBetweenness Scores");
    bc_scores = bc.scores();
    err1 = 0;
    for(count i=0; i<n; i++) {
        double x = dynbc_scores[i]-bc_scores[i];
        if (x > err1)
            err1 = x;
    }
    DEBUG("After the edge insertion: ");
}


TEST_F(DynBetweennessGTest, runApproxBetweenness) {
    //METISGraphReader reader;
    DorogovtsevMendesGenerator generator(100);
    Graph G1 = generator.generate();
    Graph G(G1, true, false);
    ApproxBetweenness bc(G, 0.1, 0.1);
    bc.run();
    DEBUG("Number of samples: ", bc.numberOfSamples());
    ApproxBetweenness bc1(G1, 0.1, 0.1);
    bc1.run();
    DEBUG("Number of samples: ", bc1.numberOfSamples());
}

TEST_F(DynBetweennessGTest, runDynVsStaticCaseInsertDirected){
    Aux::Random::setSeed(0, false);

    auto genEdgeInsert = [](const Graph& g){
        node u = GraphTools::randomNode(g);
        node v = GraphTools::randomNode(g);
        while (u == v || g.hasEdge(u, v)){
            u = GraphTools::randomNode(g);
            v = GraphTools::randomNode(g);
        }
        return std::make_pair(u, v);
    };

    for(count n = 2; n <= 25; n++)
        for(count t = 0; t < 100; t++){
            auto g = ErdosRenyiGenerator(n, 0.3, true).generate();

            if(g.numberOfEdges() == g.numberOfNodes() * (g.numberOfNodes() - 1))continue;

            auto e = genEdgeInsert(g);
            auto ibet = DynBetweenness(g);

            ibet.run();

            auto ne = GraphEvent(GraphEvent::EDGE_ADDITION, e.first, e.second);
            g.addEdge(e.first, e.second);
            ibet.update(ne);

            EXPECT_TRUE(g.hasEdge(e.first, e.second));
            auto brandes = Betweenness(g);
            brandes.run();
            g.forNodes([&](node & v){
                EXPECT_NEAR(brandes.score(v), ibet.score(v), 1e-8);
            });
        }
}

TEST_F(DynBetweennessGTest, runDynVsStaticCaseInsertUndirected){
    Aux::Random::setSeed(0, false);

    auto genEdgeInsert = [](const Graph& g){
        node u = GraphTools::randomNode(g);
        node v = GraphTools::randomNode(g);
        while (u == v || g.hasEdge(u, v)){
            u = GraphTools::randomNode(g);
            v = GraphTools::randomNode(g);
        }
        return std::make_pair(u, v);
    };

    for(count n = 2; n <= 25; n++)
        for(count t = 0; t < 100; t++){
            auto g = ErdosRenyiGenerator(n, 0.3, false).generate();
            if(g.numberOfEdges() == g.numberOfNodes() * (g.numberOfNodes() - 1) / 2)continue;

            auto e = genEdgeInsert(g);
            auto ibet = DynBetweenness(g);

            ibet.run();

            auto ne = GraphEvent(GraphEvent::EDGE_ADDITION, e.first, e.second);
            g.addEdge(e.first, e.second);
            ibet.update(ne);
            EXPECT_TRUE(g.hasEdge(e.first, e.second));

            auto brandes = Betweenness(g);
            brandes.run();
            g.forNodes([&](node & v){
                EXPECT_NEAR(brandes.score(v), ibet.score(v), 1e-8);
            });
        }
}


} /* namespace NetworKit */
