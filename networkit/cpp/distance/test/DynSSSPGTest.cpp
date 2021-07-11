// no-networkit-format
/*
 * dynSSSPGTest.cpp
 *
 *  Created on: 21.07.2014
 *      Author: ebergamini
 */

#include <gtest/gtest.h>

#include <networkit/distance/DynBFS.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/DynDijkstra.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <random>


namespace NetworKit {

class DynSSSPGTest: public testing::Test{};

TEST_F(DynSSSPGTest, testDynamicBFS_1edge) {
/* Graph:
    0    3   6
     \  / \ /
      2    5
     /  \ / \
    1    4   7
 */
    count n = 8;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(5, 7);

    BFS bfs(G, 0);
    bfs.run();
    DynBFS dbfs(G, 0);
    dbfs.run();
    std::vector<GraphEvent> batch(1);
    batch[0].type = GraphEvent::EDGE_ADDITION;
    batch[0].u = 0;
    batch[0].v = 6;
    batch[0].w = 1.0;
    for (GraphEvent edge : batch) {
        G.addEdge(edge.u, edge.v, edge.w);
    }
    dbfs.updateBatch(batch);
    bfs.run();
    G.forNodes([&] (node i) {
        EXPECT_EQ(bfs.distance(i), dbfs.distance(i));
        EXPECT_EQ(bfs.numberOfPaths(i), dbfs.numberOfPaths(i));
    });
}

TEST_F(DynSSSPGTest, testDynamicBFS_batch) {
/* Graph:
        0    3   6
        \  / \ /
            2    5
        /  \ / \
        1    4   7
*/
    count n = 8;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(5, 7);

    BFS bfs(G, 0);
    bfs.run();
    DynBFS dbfs(G, 0);
    dbfs.run();
    std::vector<GraphEvent> batch(3);
    batch[0].type = GraphEvent::EDGE_ADDITION;
    batch[0].u = 3;
    batch[0].v = 7;
    batch[0].w = 1.0;
    batch[1].type = GraphEvent::EDGE_ADDITION;
    batch[1].u = 0;
    batch[1].v = 5;
    batch[1].w = 1.0;
    batch[2].type = GraphEvent::EDGE_ADDITION;
    batch[2].u = 2;
    batch[2].v = 7;
    batch[2].w = 1.0;
    for (GraphEvent edge : batch) {
        G.addEdge(edge.u, edge.v, edge.w);
    }
    dbfs.updateBatch(batch);
    bfs.run();
    G.forNodes([&] (node i) {
        EXPECT_EQ(bfs.distance(i), dbfs.distance(i));
        EXPECT_EQ(bfs.numberOfPaths(i), dbfs.numberOfPaths(i));
    });

}


TEST_F(DynSSSPGTest, testDynamicDijkstra) {
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

    Dijkstra dij(G, 0);
    dij.run();
    DynDijkstra ddij(G, 0);
    ddij.run();
    std::vector<GraphEvent> batch(3);
    batch[0].type = GraphEvent::EDGE_ADDITION;
    batch[0].u = 0;
    batch[0].v = 4;
    batch[0].w = 1.0;
    batch[1].type = GraphEvent::EDGE_ADDITION;
    batch[1].u = 1;
    batch[1].v = 4;
    batch[1].w = 1.0;
    batch[2].type = GraphEvent::EDGE_ADDITION;
    batch[2].u = 6;
    batch[2].v = 7;
    batch[2].w = 3.0;
    for (GraphEvent edge : batch) {
        G.addEdge(edge.u, edge.v, edge.w);
    }
    ddij.updateBatch(batch);
    dij.run();
    G.forNodes([&] (node i) {
        EXPECT_EQ(dij.distance(i), ddij.distance(i));
        EXPECT_EQ(dij.numberOfPaths(i), ddij.numberOfPaths(i));
    });

}

TEST_F(DynSSSPGTest, testDynamicBFSGeneratedGraph) {
    METISGraphReader reader;
    DorogovtsevMendesGenerator generator(500);
    Graph G = generator.generate();
    DEBUG("Generated graph of dimension ", G.upperNodeIdBound());
    DynBFS dyn_bfs(G, 0);
    BFS bfs(G, 0);
    dyn_bfs.run();
    bfs.run();
    DEBUG("Before the edge insertion: ");
    count nInsertions = 750, i = 0;
    while (i < nInsertions) {
        DEBUG("Sampling a new edge");
        node v1 = GraphTools::randomNode(G);
        node v2 = GraphTools::randomNode(G);
        if (v1 != v2 && !G.hasEdge(v1, v2)) {
            i++;
            DEBUG("Adding edge number ", i);
            G.addEdge(v1, v2);
            std::vector<GraphEvent> batch;
            batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0));
            DEBUG("Running update with dynamic bfs");
            dyn_bfs.updateBatch(batch);
            DEBUG("Running from scratch with bfs");
            bfs.run();
            G.forNodes([&] (node i) {
            //	std::cout<<"Node "<<i<<":"<<std::endl;
            //	std::cout<<"Actual distance: "<<dij.distance(i)<<", computed distance: "<<ddij.distance(i)<<std::endl;
            //	std::cout<<"Actual number of paths: "<<dij.numberOfPaths(i)<<", computed one: "<<ddij.numberOfPaths(i)<<std::endl;
                EXPECT_EQ(dyn_bfs.distance(i), bfs.distance(i));
                EXPECT_EQ(dyn_bfs.numberOfPaths(i), bfs.numberOfPaths(i));
            });
        }
    }
}

TEST_F(DynSSSPGTest, testDynamicDijkstraGeneratedGraph) {
    METISGraphReader reader;
    DorogovtsevMendesGenerator generator(1000);
    Graph G1 = generator.generate();
    Graph G(G1, true, false);
    DEBUG("Generated graph of dimension ", G.upperNodeIdBound());
    DynDijkstra dyn_dij(G, 0);
    Dijkstra dij(G, 0);
    dyn_dij.run();
    dij.run();
    DEBUG("Before the edge insertion: ");
    count nInsertions = 10, i = 0;
    while (i < nInsertions) {
        DEBUG("Sampling a new edge");
        node v1 = GraphTools::randomNode(G);
        node v2 = GraphTools::randomNode(G);
        if (v1 != v2 && !G.hasEdge(v1, v2)) {
            i++;
            DEBUG("Adding edge number ", i);
            G.addEdge(v1, v2);
            std::vector<GraphEvent> batch;
            batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, 1.0));
            DEBUG("Running update with dynamic dijkstra");
            dyn_dij.updateBatch(batch);
            DEBUG("Running from scratch with dijkstra");
            dij.run();
            G.forNodes([&] (node i) {
            //	std::cout<<"Node "<<i<<":"<<std::endl;
            //	std::cout<<"Actual distance: "<<dij.distance(i)<<", computed distance: "<<ddij.distance(i)<<std::endl;
            //	std::cout<<"Actual number of paths: "<<dij.numberOfPaths(i)<<", computed one: "<<ddij.numberOfPaths(i)<<std::endl;
                EXPECT_EQ(dyn_dij.distance(i), dij.distance(i));
                EXPECT_EQ(dyn_dij.numberOfPaths(i), dij.numberOfPaths(i));
            });
        }
    }
}

TEST_F(DynSSSPGTest, testDynamicDijkstraBatches) {
    METISGraphReader reader;
    std::default_random_engine random_generator;
    std::normal_distribution<double> distribution(100,10);
    DorogovtsevMendesGenerator generator(100);
    Graph G1 = generator.generate();
    Graph G(G1, true, false);
    DEBUG("Generated graph of dimension ", G.upperNodeIdBound());
    // add random normal weights to G

    G.forNodes([&] (node source) {
        DynDijkstra dyn_dij(G, source, true);
        Dijkstra dij(G, source);
        dyn_dij.run();
        dij.run();
        DEBUG("Before the edge insertion: ");
        count batchSize = 8;
        count nBatches = 1, i = 0;
        for (count j=0; j<nBatches; j++) {
            std::vector<GraphEvent> batch;
            i = 0;
            while (i < batchSize) {
                DEBUG("Sampling a new edge");
                node v1 = GraphTools::randomNode(G);
                node v2 = GraphTools::randomNode(G);
                if (v1 != v2 && !G.hasEdge(v1, v2)) {
                    i++;
                    double number = distribution(random_generator);
                    G.addEdge(v1, v2, number);
                    batch.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, v1, v2, number));
                }
            }
            DEBUG("batch size: ", batch.size());
            DEBUG("Updating with dynamic dijkstra");
            dyn_dij.updateBatch(batch);
            DEBUG("Running from scratch with dijkstra");
            dij.run();
            G.forNodes([&] (node i) {
            //	std::cout<<"Node "<<i<<":"<<std::endl;
            //	std::cout<<"Actual distance: "<<dij.distance(i)<<", computed distance: "<<ddij.distance(i)<<std::endl;
            //	std::cout<<"Actual number of paths: "<<dij.numberOfPaths(i)<<", computed one: "<<ddij.numberOfPaths(i)<<std::endl;
                EXPECT_EQ(dyn_dij.distance(i), dij.distance(i));
                EXPECT_EQ(dyn_dij.numberOfPaths(i), dij.numberOfPaths(i));
                if (i != source)
                    assert(dyn_dij.distance(i) != 0);
            //	EXPECT_EQ(dyn_dij.getPredecessors(i).size(), dij.getPredecessors(i).size());
            });
        }
    });
}

} /* namespace NetworKit */
