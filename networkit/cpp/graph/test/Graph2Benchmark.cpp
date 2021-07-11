// no-networkit-format
/*
 * Graph2Benchmark.cpp
 *
 *  Created on: 05.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/auxiliary/Log.hpp>

namespace NetworKit {

class Graph2Benchmark: public testing::Test {};

TEST_F(Graph2Benchmark, graphConstruction) {
    count n = 1e+7;;

    Aux::Timer run;
    INFO("[BEGIN] (n=" , n , ")");
    run.start();

    Graph G(n);

    run.stop();
    INFO("[DONE]" , run.elapsedTag());
}


TEST_F(Graph2Benchmark, nodeIteration) {
    count n = 1e+7;;
    Graph G(n);


    std::vector<node> nodes(n, 0);

    Aux::Timer run;
    INFO("[BEGIN] (n=" , n , ")");
    run.start();

    G.forNodes([&](node v){
        nodes[v] = v;
    });

    run.stop();
    INFO("[DONE]" , run.elapsedTag());

}

TEST_F(Graph2Benchmark, parallelNodeIteration) {
    count n = 1e+7;;
    Graph G(n);


    std::vector<node> nodes(n, 0);

    Aux::Timer run;
    INFO("[BEGIN] (n=" , n , ")");
    run.start();

    G.parallelForNodes([&](node v){
        nodes[v] = v;
    });

    run.stop();
    INFO("[DONE]" , run.elapsedTag());
}


TEST_F(Graph2Benchmark, nodePairIteration) {
    count n = 1e+4;;
    Graph G(n);


    Aux::Timer run;
    INFO("[BEGIN] (n=" , n , ")");
    run.start();

    count p = 0;
    G.forNodePairs([&](node, node){
        p += 1;
    });

    run.stop();
    INFO("[DONE]" , run.elapsedTag());

    EXPECT_EQ((n * (n-1)) / 2, p);
}


TEST_F(Graph2Benchmark, edgeInsertion) {

    count n = 1e+4;
    Graph G(n);

    Aux::Timer run;
    INFO("[BEGIN] (n=" , n , ")");
    run.start();

    G.forNodePairs([&](node u, node v){
        G.addEdge(u, v);
    });


    run.stop();
    INFO("[DONE]" , run.elapsedTag());

}


TEST_F(Graph2Benchmark, edgeRemoval) {
    count n = 1e+4;
    Graph G(n);

    // insert edges
    G.forNodePairs([&](node u, node v){
        G.addEdge(u, v);
    });


    Aux::Timer run;
    INFO("[BEGIN] (n=" , n , ")");
    run.start();

    G.forNodePairs([&](node u, node v){
        G.removeEdge(u, v);
    });

    run.stop();
    INFO("[DONE]" , run.elapsedTag());

}


TEST_F(Graph2Benchmark, edgeIteration) {
    count n = 1e+4;
    Graph G(n);

    // insert edges
    G.forNodePairs([&](node u, node v){
        G.addEdge(u, v);
    });


    Aux::Timer run;
    INFO("[BEGIN] (n=" , n , ")");
    run.start();

    G.forEdges([&](node, node){

    });

    run.stop();
    INFO("[DONE]" , run.elapsedTag());
}

TEST_F(Graph2Benchmark, parallelEdgeIteration) {
    count n = 1e+4;
    Graph G(n);

    // insert edges
    G.forNodePairs([&](node u, node v){
        G.addEdge(u, v);
    });

    Aux::Timer run;
    INFO("[BEGIN] (n=" , n , ")");
    run.start();

    count i = 0;
    G.parallelForEdges([&](node, node){
        i += 1;
    });

    EXPECT_TRUE(true) << "just iterate";

    run.stop();
    INFO("[DONE]" , run.elapsedTag());
}

TEST_F(Graph2Benchmark, parallelSumForNodes) {
    count n = 1e+7;
    Graph G(n);

    double sum = G.parallelSumForNodes([&](node) {
        return 1;
    });

    EXPECT_EQ(n, sum) << "each node adds 1";

}



TEST_F(Graph2Benchmark, nodeInsertion) {
    count n = 1e+4;

    Graph G(0); // empty graph

    for (count i = 0; i < n; ++i) {
        //node v = G.addNode();
        G.addNode();
    }

    EXPECT_EQ(n, G.numberOfNodes()) << "n nodes should have been added";
}

TEST_F(Graph2Benchmark, nodeRemoval) {
    count n = 1e+4;

    Graph G(n); // empty graph

    for (node u = 0; u < n; ++u) {
        G.removeNode(u);
    }

    EXPECT_EQ(0u, G.numberOfNodes()) << "no nodes should be left";
}

} /* namespace NetworKit */
