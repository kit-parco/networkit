// no-networkit-format
/*
 * APSPGTest.cpp
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe, Elisabetta Bergamini
 */

#include <string>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/distance/APSP.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/DynAPSP.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {
class APSPGTest: public testing::Test{};

TEST_F(APSPGTest, testAPSP) {
/* Graph:
     ______
        /      \
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
    G.addEdge(0, 6);

    APSP apsp(G);
    apsp.run();
    std::vector<std::vector<edgeweight> > distances = apsp.getDistances();
    ASSERT_EQ(distances[0][0],0);
    ASSERT_EQ(distances[0][1],2);
    ASSERT_EQ(distances[0][2],1);
    ASSERT_EQ(distances[0][3],2);
    ASSERT_EQ(distances[0][4],2);
    ASSERT_EQ(distances[0][5],2);
    ASSERT_EQ(distances[0][6],1);
    ASSERT_EQ(distances[1][0],2);
    ASSERT_EQ(distances[1][1],0);
    ASSERT_EQ(distances[1][2],1);
    ASSERT_EQ(distances[1][3],2);
    ASSERT_EQ(distances[1][4],2);
    ASSERT_EQ(distances[1][5],3);
    ASSERT_EQ(distances[1][6],3);
}

TEST_F(APSPGTest, testAPSPUnweightedER) {
    Aux::Random::setSeed(42, false);
    auto G = ErdosRenyiGenerator(100, 0.01, false).generate();
    G.addNode();

    APSP apsp(G);
    apsp.run();
    G.forNodes([&](const node u) {
        BFS bfs(G, u, false);
        bfs.run();
        const auto &dist = bfs.getDistances();
        G.forNodes([&](const node v) {
            EXPECT_DOUBLE_EQ(dist[v], apsp.getDistance(u, v));
        });
    });
}

TEST_F(APSPGTest, testAPSPWeightedER) {
    Aux::Random::setSeed(42, false);
    auto G = ErdosRenyiGenerator(100, 0.01, false).generate();
    G.addNode();
    G = GraphTools::toWeighted(G);

    APSP apsp(G);
    apsp.run();
    G.forNodes([&](const node u) {
        Dijkstra dij(G, u, false);
        dij.run();
        const auto &dist = dij.getDistances();
        G.forNodes([&](const node v) {
            EXPECT_DOUBLE_EQ(dist[v], apsp.getDistance(u, v));
        });
    });
}

TEST_F(APSPGTest, debugAPSP) {
    count n = 1000;
    count m = int(n * n);
    Graph G(n, true, false);
    for (count i = 0; i < m; i++) {
        node u = GraphTools::randomNode(G);
        node v = GraphTools::randomNode(G);
        if (u != v && !G.hasEdge(u, v)) {
            G.addEdge(u, v, Aux::Random::integer(10));
        }
    }
    INFO("Nodes: ", G.numberOfNodes(), ", edges: ", G.numberOfEdges());
    APSP apsp(G);
    apsp.run();
}

TEST_F(APSPGTest, testDynAPSPRealGraph) {
    METISGraphReader reader;
    Graph G = reader.read("input/karate.graph");
    DynAPSP apsp(G);
    apsp.run();
    for (count i = 0; i < 10; i++) {
        count u = GraphTools::randomNode(G);
        count v = GraphTools::randomNode(G);
        while(G.hasEdge(u, v)) {
            u = GraphTools::randomNode(G);
            v = GraphTools::randomNode(G);
        }
        DEBUG("u = ", u, ", v = ", v);
        GraphEvent event(GraphEvent::EDGE_ADDITION, u, v, 1);
        G.addEdge(u, v, 1);
        apsp.update(event);
        DynAPSP apsp2(G);
        apsp2.run();
        std::vector<std::vector<edgeweight> > distances = apsp.getDistances();
        std::vector<std::vector<edgeweight> > distances2 = apsp2.getDistances();
        G.forNodes([&](node i) {
            G.forNodes([&](node j) {
                if (distances[i][j] != distances2[i][j])
                    DEBUG("i, j = ", i, " ", j, ", dist[i][j] = ", distances[i][j], ", dist2[i][j] = ", distances2[i][j], ", u = ", u, ", v = ", v);
                EXPECT_NEAR(distances[i][j], distances2[i][j], 0.0001);
            });
        });
 }
}

TEST_F(APSPGTest, testAPSPRunDirectedUnweighted) {
    // build G
    Graph G(6, false, true);
    G.addEdge(0,1);
    G.addEdge(1,2);
    G.addEdge(2,3);
    G.addEdge(3,4);
    G.addEdge(4,5);
    G.addEdge(5,0);
    DynAPSP apsp(G);
    apsp.run();
    std::vector<std::vector<edgeweight> > distances = apsp.getDistances();

    EXPECT_EQ(distances[0][0], 0);
    EXPECT_EQ(distances[0][2], 2);
    EXPECT_EQ(distances[0][5], 5);

}

TEST_F(APSPGTest, testAPSPRunUndirectedWeighted) {
    // build G
    Graph G(4, true, false); // set to 5 and test for isolated nodes and disconnected components
    G.addEdge(0,1,1);
    G.addEdge(0,2,3);
    G.addEdge(1,2,1);
    G.addEdge(1,3,1);
    DynAPSP apsp(G);
    apsp.run();
    std::vector<std::vector<edgeweight> > distances = apsp.getDistances();

    EXPECT_EQ(distances[0][0], 0);
    EXPECT_EQ(distances[0][2], 2);
    EXPECT_EQ(distances[0][3], 2);
}

TEST_F(APSPGTest, testAPSPInsertionUndirectedUnweighted) {
    // build G
    // 0 --- 1 --- 2
  //       |
  //       3 --- 4 --- 5 --- 6

    Graph G(7);
    G.addEdge(0,1);
    G.addEdge(1,2);
    G.addEdge(1,3);
    G.addEdge(3,4);
    G.addEdge(4,5);
    G.addEdge(5,6);

    // run dyn apsp with insertion of edge (2, 5)
    DynAPSP apsp(G);
    apsp.run();
    std::vector<std::vector<edgeweight> > distances = apsp.getDistances();

    // apply graph update
    G.addEdge(2, 5);
    GraphEvent event(GraphEvent::EDGE_ADDITION, 2, 5, 1);
    apsp.update(event);
    distances = apsp.getDistances();

    apsp.run();
    std::vector<std::vector<edgeweight> > distances2 = apsp.getDistances();

    G.forNodes([&](node i) {
        G.forNodes([&](node j) {
            if (distances[i][j] != distances2[i][j]){
                DEBUG("i = ", i, ", j = ", j, ", distance: ", distances[i][j], ", expected: ", distances2[i][j]);
            }
            EXPECT_EQ(distances[i][j], distances2[i][j]);
        });
    });
}

TEST_F(APSPGTest, testAPSPInsertionDirectedUnweighted) {
    Graph G(6, false, true);
    G.addEdge(0,1);
    G.addEdge(1,2);
    G.addEdge(2,3);
    G.addEdge(3,4);
    G.addEdge(4,5);
    G.addEdge(5,0);
    DynAPSP apsp(G);
    apsp.run();
    std::vector<std::vector<edgeweight> > distances = apsp.getDistances();

    G.addEdge(0,3);
    GraphEvent event(GraphEvent::EDGE_ADDITION, 0, 3, 1);
    apsp.update(event);

    distances = apsp.getDistances();

    DynAPSP APSP(G);
    APSP.run();
    std::vector<std::vector<edgeweight> > distances2 = APSP.getDistances();

    G.forNodes([&](node i) {
        G.forNodes([&](node j) {
            EXPECT_EQ(distances[i][j], distances2[i][j]) << "i, j = " << i << ", " << j;
        });
    });

    std::vector<node> path10 = apsp.getPath(1, 0);
    DEBUG("path10 = ", path10);
}

TEST_F(APSPGTest, testAPSPUndirectedWeighted) {
    // build G
    Graph G(7, true, false);
    G.addEdge(0,1, 1);
    G.addEdge(1,2, 0.01);
    G.addEdge(1,3, 0.1);
    G.addEdge(3,4, 0.001);
    G.addEdge(4,5, 0.0001);
    G.addEdge(5,6, 0.00001);

    // Run baseline apsp with ID 1
    DynAPSP apsp(G);
    apsp.run();
    std::vector<std::vector<edgeweight> > distances = apsp.getDistances();

    // apply graph update edge addition with ID 2
    G.addEdge(2, 5, 0.000002);
    GraphEvent event2(GraphEvent::EDGE_ADDITION, 2, 5, 0.000002);
    apsp.update(event2);

    distances = apsp.getDistances();

    DynAPSP APSP(G);
    APSP.run();
    std::vector<std::vector<edgeweight> > distances2 = APSP.getDistances();

    G.forNodes([&](node i) {
        G.forNodes([&](node j) {
            EXPECT_NEAR(distances[i][j], distances2[i][j], 0.0001) << "i, j = " << i << ", " << j;
        });
    });

    // apply graph update edge weight update with ID 4
    G.setWeight(2, 5, 0.000001);
    GraphEvent event4(GraphEvent::EDGE_WEIGHT_INCREMENT, 2, 5, -0.000001);
    apsp.update(event4);
    distances = apsp.getDistances();

    DynAPSP apsp4(G);
    apsp4.run();
    std::vector<std::vector<edgeweight> > distances4 = apsp4.getDistances();

    G.forNodes([&](node i) {
        G.forNodes([&](node j) {
            EXPECT_NEAR(distances[i][j], distances4[i][j], 0.0001) << "i, j = " << i << ", " << j;
        });
    });
}

TEST_F(APSPGTest, testAPSPDirectedWeighted1) {
    Graph G(5, true, true); // G+ Ghouse
    G.addEdge(3,1,1);
    G.addEdge(1,0,2);
    G.addEdge(0,2,3);
    G.addEdge(2,1,4);
    G.addEdge(1,4,5);
    G.addEdge(4,3,6);
    G.addEdge(3,2,7);
    G.addEdge(2,4,8);

    // Run baseline apsp with ID 1
    DynAPSP apsp(G);
    apsp.run();
    std::vector<std::vector<edgeweight> > distances = apsp.getDistances();

    // apply graph update edge insertion update with ID 2
    //INFO("entering update 2");
    G.addEdge(3, 1, 1);
    GraphEvent event2(GraphEvent::EDGE_ADDITION, 3, 1, 1);
    apsp.update(event2);

    distances = apsp.getDistances();

    DynAPSP APSP(G);
    APSP.run();
    std::vector<std::vector<edgeweight> > distances2 = APSP.getDistances();
    G.forNodes([&](node i) {
        G.forNodes([&](node j) {
            EXPECT_NEAR(distances[i][j], distances2[i][j], 0.0001) << "i, j = " << i << ", " << j;
        });
    });
}

TEST_F(APSPGTest, testAPSPDirectedWeighted2) {
    Graph G(6, true, true);
    G.addEdge(0,1,1);
    G.addEdge(1,2,1);
    G.addEdge(2,3,1);
    G.addEdge(3,4,1);
    G.addEdge(4,5,1);
    G.addEdge(3,5,3);

    // Run baseline apsp with ID 1
    DynAPSP apsp(G);
    apsp.run();
    std::vector<std::vector<edgeweight> > distances = apsp.getDistances();

    // apply graph update edge insertion update with ID 2
    //INFO("entering update 2");
    G.addEdge(1, 3, 1);
    GraphEvent event2(GraphEvent::EDGE_ADDITION, 1, 3, 1);
    apsp.update(event2);

    distances = apsp.getDistances();

    DynAPSP APSP(G);
    APSP.run();
    std::vector<std::vector<edgeweight> > distances2 = APSP.getDistances();
    G.forNodes([&](node i) {
        G.forNodes([&](node j) {
            EXPECT_NEAR(distances[i][j], distances2[i][j], 0.0001) << "i, j = " << i << ", " << j;
        });
    });
}

TEST_F(APSPGTest, testAPSPIsolatedNode) {
    Graph G(3, false, false);

    // Run baseline apsp with ID 1
    DynAPSP apsp(G);
    apsp.run();
    std::vector<std::vector<edgeweight> > distances = apsp.getDistances();

    // apply graph update edge insertion update with ID 2
    //INFO("entering update 2");
    G.addEdge(0,1);
    GraphEvent event2(GraphEvent::EDGE_ADDITION, 0, 1, 1);
    apsp.update(event2);

    distances = apsp.getDistances();

    DynAPSP APSP(G);
    APSP.run();
    std::vector<std::vector<edgeweight> > distances2 = APSP.getDistances();
    G.forNodes([&](node i) {
        G.forNodes([&](node j) {
            EXPECT_NEAR(distances[i][j], distances2[i][j], 0.0001);
        });
    });

    std::vector<node> path;
    EXPECT_EQ(apsp.getPath(0, 2), path);
}

TEST_F(APSPGTest, testAPSPEventTypeError) {
    Graph G(6, true, true);
    G.addEdge(0,1,1);
    G.addEdge(1,2,1);
    G.addEdge(2,3,1);
    G.addEdge(3,4,1);
    G.addEdge(4,5,1);
    G.addEdge(3,5,3);

    // Run baseline apsp with ID 1
    DynAPSP apsp(G);
    apsp.run();
    std::vector<std::vector<edgeweight> > distances = apsp.getDistances();
    //INFO("distances[0]: ", distances[0][0], distances[0][1], distances[0][2], distances[0][3], distances[0][4], distances[0][5], distances[0][6]);

    // test throw 1. Edge deletions are not allowed.
    G.removeEdge(0, 1);
    GraphEvent event2(GraphEvent::EDGE_REMOVAL, 0, 1, 1);
    EXPECT_ANY_THROW(apsp.update(event2));

    // test throw 2. Edge weight increases are not allowed
    G.setWeight(1, 2, 2);
    GraphEvent event3(GraphEvent::EDGE_WEIGHT_INCREMENT, 1, 2, 1);
    EXPECT_ANY_THROW(apsp.update(event3));
}



} /* namespace NetworKit */
