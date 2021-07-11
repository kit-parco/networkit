// no-networkit-format
/*
* ConnectedComponentsGTest.cpp
*
*  Created on: Sep 16, 2013
*      Author: Maximilian Vogel
*/
#include <gtest/gtest.h>

#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/components/ParallelConnectedComponents.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/components/DynConnectedComponents.hpp>
#include <networkit/components/DynWeaklyConnectedComponents.hpp>
#include <networkit/components/WeaklyConnectedComponents.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/BFS.hpp>

#include <networkit/distance/Diameter.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/io/KONECTGraphReader.hpp>
#include <networkit/generators/HavelHakimiGenerator.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>

namespace NetworKit {

class ConnectedComponentsGTest: public testing::Test{};

TEST_F(ConnectedComponentsGTest, testConnectedComponentsTiny) {
    // construct graph
    Graph g(20);

    g.addEdge(0,1,0);
    g.addEdge(1,2,0);
    g.addEdge(2,4,0);
    g.addEdge(4,8,0);
    g.addEdge(8,16,0);
    g.addEdge(16,19,0);

    g.addEdge(3,5,0);
    g.addEdge(5,6,0);
    g.addEdge(6,7,0);
    g.addEdge(7,9,0);

    g.addEdge(10,11,0);
    g.addEdge(10,18,0);
    g.addEdge(10,12,0);
    g.addEdge(18,17,0);

    g.addEdge(13,14,0);

    // initialize ConnectedComponents
    ConnectedComponents ccs(g);
    ccs.run();

    // check result
    EXPECT_EQ(5, ccs.numberOfComponents());
    EXPECT_TRUE(ccs.componentOfNode(0) == ccs.componentOfNode(19));
    EXPECT_TRUE(ccs.componentOfNode(3) == ccs.componentOfNode(7));
}


TEST_F(ConnectedComponentsGTest, testConnectedComponents) {
    // construct graph
    METISGraphReader reader;
    Graph G = reader.read("input/astro-ph.graph");
    ConnectedComponents cc(G);
    cc.run();
    DEBUG("Number of components: ", cc.numberOfComponents());
    EXPECT_EQ(1029u, cc.numberOfComponents());
}

TEST_F(ConnectedComponentsGTest, testParallelConnectedComponents) {
    METISGraphReader reader;
    std::vector<std::string> graphs = {"PGPgiantcompo", "celegans_metabolic", "hep-th", "jazz"};

    for (auto graphName: graphs) {
        Graph G = reader.read("input/" + graphName + ".graph");
        ParallelConnectedComponents cc(G);
        cc.runSequential();
        count seqNum = cc.numberOfComponents();
        cc.run();
        count parNum = cc.numberOfComponents();
        DEBUG("Number of components: ", seqNum);
        EXPECT_EQ(seqNum, parNum);
    }
}

TEST_F(ConnectedComponentsGTest, testParallelConnectedComponentsWithDeletedNodes) {
    Graph G(100);
    G.forNodePairs([&](node u, node v){
        G.addEdge(u,v);
    });


    {
        ParallelConnectedComponents cc(G);
        cc.run();
        EXPECT_EQ(1, cc.numberOfComponents()) << "The complete graph has just one connected component";
    }

    for (node u = 0; u < 10; ++u) {
        G.forNeighborsOf(u, [&](node v) {
            G.removeEdge(u, v);
        });
        G.removeNode(u);
    }

    {
        ParallelConnectedComponents cc(G);
        cc.run();
        EXPECT_EQ(1, cc.numberOfComponents()) << "The complete graph with 10 nodes removed has still just one connected component (run())";
        cc.runSequential();
        EXPECT_EQ(1, cc.numberOfComponents()) << "The complete graph with 10 nodes removed has still just one connected component (runSequential())";
    }

}

TEST_F(ConnectedComponentsGTest, benchConnectedComponents) {
    // construct graph
    METISGraphReader reader;
    Graph G = reader.read("input/coAuthorsDBLP.graph");
    ConnectedComponents cc(G);
    cc.run();
    DEBUG("Number of components: ", cc.numberOfComponents());
    EXPECT_EQ(1u, cc.numberOfComponents());
}


TEST_F(ConnectedComponentsGTest, testStronglyConnectedComponentsTiny) {

    auto comparePartitions = [](const Partition& p1, const Partition& p2) {
        std::vector<index> partitionIdMap(p1.upperBound(), none);
        ASSERT_EQ(p1.numberOfElements(), p2.numberOfElements());
        ASSERT_EQ(p1.numberOfSubsets(), p2.numberOfSubsets());

        p1.forEntries([&](node v, index p) {
            if (partitionIdMap[p] == none) {
                partitionIdMap[p] = p2.subsetOf(v);
            }
            index p_mapped = partitionIdMap[p];
            ASSERT_EQ(p_mapped, p);
        });
    };


    count n = 8;
    count m = 14;
    Graph G(n, false, true);

    G.addEdge(0, 4);
    G.addEdge(1, 0);
    G.addEdge(2, 1);
    G.addEdge(2, 3);
    G.addEdge(3, 2);
    G.addEdge(4, 1);
    G.addEdge(5, 1);
    G.addEdge(5, 4);
    G.addEdge(5, 6);
    G.addEdge(6, 2);
    G.addEdge(6, 5);
    G.addEdge(7, 3);
    G.addEdge(7, 6);
    G.addEdge(7, 7);

    ASSERT_EQ(n, G.numberOfNodes());
    ASSERT_EQ(m, G.numberOfEdges());

    count z = G.upperNodeIdBound();
    Partition p_expected(z);
    p_expected.allToSingletons();
    p_expected[0] = 0;
    p_expected[1] = 0;
    p_expected[2] = 1;
    p_expected[3] = 1;
    p_expected[4] = 0;
    p_expected[5] = 2;
    p_expected[6] = 2;
    p_expected[7] = 3;
    p_expected.compact();

    StronglyConnectedComponents scc(G);
    scc.run();
    Partition p_actual = scc.getPartition();
    p_actual.compact();

    comparePartitions(p_expected, p_actual);
}

TEST_F(ConnectedComponentsGTest, testStronglyConnectedComponents) {

    auto testComponent = [](const Graph &G, const std::vector<node> &cmp) {
        std::vector<bool> inComponent(G.upperNodeIdBound());
        std::vector<bool> reachableFromComponent(G.upperNodeIdBound());
        for (const auto u : cmp) {
            inComponent[u] = true;
            std::unordered_set<node> unvisited(cmp.begin(), cmp.end());
            Traversal::BFSfrom(G, u, [&](const node v, count) {
                const auto iter = unvisited.find(v);
                if (iter != unvisited.end()) {
                    unvisited.erase(iter);
                }
                reachableFromComponent[v] = true;
            });

            EXPECT_TRUE(unvisited.empty());
        }

        G.forNodes([&](const node u) {

            if (inComponent[u] || !reachableFromComponent[u])
                return;

            Traversal::BFSfrom(G, u, [&](const node v) {
                EXPECT_FALSE(inComponent[v]);
            });
        });
    };

    constexpr count n = 200;
    for (int i : {1, 2, 3}) {
        Aux::Random::setSeed(i, false);
        for (double p : {0.01, 0.05, 0.1}) {
            const auto G = ErdosRenyiGenerator(n, p, true).generate();
            StronglyConnectedComponents scc(G);
            scc.run();

            const auto nComponents = scc.numberOfComponents();
            const auto cmpVec = scc.getComponents();

            EXPECT_EQ(cmpVec.size(), nComponents);

            const auto compSizes = scc.getComponentSizes();
            for (const auto &entry : compSizes) {
                EXPECT_TRUE(entry.first < cmpVec.size());
                EXPECT_EQ(cmpVec[entry.first].size(), entry.second);
            }

            for (index i = 0; i < cmpVec.size(); ++i) {
                for (const auto u : cmpVec[i]) {
                    EXPECT_EQ(scc.componentOfNode(u), i);
                }
            }

            for (const auto &cmp : cmpVec) {
                testComponent(G, cmp);
            }
        }
    }
}

TEST_F(ConnectedComponentsGTest, testDynConnectedComponentsTiny) {
    // construct graph
    Graph g;
    for (count i = 0; i < 20; i++) {
        g.addNode();
    }
    g.addEdge(0,1,0);
    g.addEdge(1,2,0);
    g.addEdge(2,4,0);
    g.addEdge(4,8,0);
    g.addEdge(8,16,0);
    g.addEdge(16,19,0);

    g.addEdge(3,5,0);
    g.addEdge(5,6,0);
    g.addEdge(6,7,0);
    g.addEdge(7,9,0);

    g.addEdge(10,11,0);
    g.addEdge(10,18,0);
    g.addEdge(10,12,0);
    g.addEdge(18,17,0);

    g.addEdge(13,14,0);

    // initialize DynConnectedComponents
    DynConnectedComponents dccs(g);
    dccs.run();

    // check result
    EXPECT_EQ(5, dccs.numberOfComponents());
    EXPECT_TRUE(dccs.componentOfNode(0) == dccs.componentOfNode(19));
    EXPECT_TRUE(dccs.componentOfNode(3) == dccs.componentOfNode(7));

    g.addEdge(13, 15, 0);
    dccs.update(GraphEvent(GraphEvent::EDGE_ADDITION, 13, 15, 0));
    EXPECT_EQ(4, dccs.numberOfComponents());
    EXPECT_TRUE(dccs.componentOfNode(14) == dccs.componentOfNode(15));
    EXPECT_TRUE(dccs.componentOfNode(15) != dccs.componentOfNode(0));

    // Create batch and update
    std::vector<GraphEvent> batch {
        GraphEvent(GraphEvent::EDGE_ADDITION, 15, 19, 0),
        GraphEvent(GraphEvent::EDGE_REMOVAL, 6, 7, 0),
        GraphEvent(GraphEvent::EDGE_ADDITION, 7, 17, 0),
        GraphEvent(GraphEvent::EDGE_REMOVAL, 2, 4, 0)
    };

    for (auto e : batch) {
        if (e.type == GraphEvent::EDGE_ADDITION){
            g.addEdge(e.u, e.v);
        }
        if (e.type == GraphEvent::EDGE_REMOVAL) {
            g.removeEdge(e.u, e.v);
        }
    }

    dccs.updateBatch(batch);
    EXPECT_EQ(4, dccs.numberOfComponents());
    EXPECT_TRUE(dccs.componentOfNode(0) != dccs.componentOfNode(14));
    EXPECT_TRUE(dccs.componentOfNode(9) == dccs.componentOfNode(11));
    EXPECT_TRUE(dccs.componentOfNode(3) != dccs.componentOfNode(9));
}


TEST_F(ConnectedComponentsGTest, testDynConnectedComponents) {
    // construct graph
    METISGraphReader reader;
    Graph G = reader.read("input/PGPgiantcompo.graph");
    DynConnectedComponents dccs(G);
    dccs.run();
    ConnectedComponents cc(G);
    cc.run();
    DEBUG("Number of components: ", dccs.numberOfComponents());
    // Testing static run.
    EXPECT_EQ(cc.numberOfComponents(), dccs.numberOfComponents());

    int numberOfTests = 500;

    // Probability to perform an edge insertion or removal.
    float p = 0.5;
    srand (time(NULL));
    for (int i = 0; i < numberOfTests; ++i) {
        node u = 0;
        node v = 1;

        // Perform edge insertion
        if (((double) rand() / (RAND_MAX)) > p) {
            while (G.hasEdge(u, v)) {
                u = GraphTools::randomNode(G);
                v = GraphTools::randomNode(G);
            }
            G.addEdge(u, v);
            dccs.update(GraphEvent(GraphEvent::EDGE_ADDITION, u, v, 0));
        }
        else {
            while (!G.hasEdge(u, v)) {
                std::pair<node, node> edge = GraphTools::randomEdge(G);
                u = edge.first;
                v = edge.second;
            }
            G.removeEdge(u, v);
            dccs.update(GraphEvent(GraphEvent::EDGE_REMOVAL, u, v, 0));
        }
    }

    cc.run();
    EXPECT_EQ(cc.numberOfComponents(), dccs.numberOfComponents());

    // Testing batch update.
    std::vector<GraphEvent> batch(numberOfTests);
    for (int i = 0; i < numberOfTests; ++i) {
        node u = GraphTools::randomNode(G);
        node v = GraphTools::randomNode(G);
        if (((double) rand() / (RAND_MAX)) > -1) {
            while (G.hasEdge(u, v)) {
                u = GraphTools::randomNode(G);
                v = GraphTools::randomNode(G);
            }
            batch[i] = GraphEvent(GraphEvent::EDGE_ADDITION, u, v, 0);
            G.addEdge(u, v);
        }
        else {
            while (!G.hasEdge(u, v)) {
                std::pair<node, node> edge = GraphTools::randomEdge(G);
                u = edge.first;
                v = edge.second;
            }
            batch[i] = GraphEvent(GraphEvent::EDGE_REMOVAL, u, v, 0);
            G.removeEdge(u, v);
        }
    }

    dccs.updateBatch(batch);
    cc.run();
    EXPECT_EQ(cc.numberOfComponents(), dccs.numberOfComponents());
}


TEST_F(ConnectedComponentsGTest, testWeaklyConnectedComponentsTiny) {
    // construct graph
    Graph g(0, false, true);
    for (count i = 0; i < 20; ++i) {
        g.addNode();
    }
    g.addEdge(0,1,0);
    g.addEdge(1,2,0);
    g.addEdge(2,4,0);
    g.addEdge(4,8,0);
    g.addEdge(8,16,0);
    g.addEdge(16,19,0);

    g.addEdge(3,5,0);
    g.addEdge(5,6,0);
    g.addEdge(6,7,0);
    g.addEdge(7,9,0);

    g.addEdge(10,11,0);
    g.addEdge(10,18,0);
    g.addEdge(10,12,0);
    g.addEdge(18,17,0);
    g.addEdge(17,10,0);

    g.addEdge(13,14,0);

    // initialize WeaklyConnectedComponents
    WeaklyConnectedComponents wcc(g);
    wcc.run();

    // check result
    EXPECT_EQ(5, wcc.numberOfComponents());
    EXPECT_TRUE(wcc.componentOfNode(0) == wcc.componentOfNode(19));
    EXPECT_TRUE(wcc.componentOfNode(3) == wcc.componentOfNode(7));
}

TEST_F(ConnectedComponentsGTest, testWeaklyConnectedComponents) {
    // construct graph
    EdgeListReader directReader(' ', 0, "%", false, true);
    Graph G = directReader.read("input/johnson8-4-4.edgelist");
    EdgeListReader undirectReader(' ', 0, "%", false, false);
    Graph Gu = undirectReader.read("input/johnson8-4-4.edgelist");
    WeaklyConnectedComponents wc(G);
    ConnectedComponents cc(Gu);
    wc.run();
    cc.run();
    DEBUG("Number of components: ", cc.numberOfComponents());
    EXPECT_EQ(wc.numberOfComponents(), cc.numberOfComponents());
}


TEST_F(ConnectedComponentsGTest, testDynWeaklyConnectedComponentsTiny) {
    // construct graph
    Graph g(20, false, true);

    g.addEdge(0,1,0);
    g.addEdge(1,2,0);
    g.addEdge(2,4,0);
    g.addEdge(4,8,0);
    g.addEdge(8,16,0);
    g.addEdge(16,19,0);

    g.addEdge(3,5,0);
    g.addEdge(5,6,0);
    g.addEdge(6,7,0);
    g.addEdge(7,9,0);

    g.addEdge(10,11,0);
    g.addEdge(10,18,0);
    g.addEdge(10,12,0);
    g.addEdge(18,17,0);

    g.addEdge(13,14,0);

    // initialize DynWeaklyConnectedComponents
    DynWeaklyConnectedComponents dw(g);
    dw.run();

    // check result
    EXPECT_EQ(5, dw.numberOfComponents());
    EXPECT_EQ(dw.componentOfNode(0), dw.componentOfNode(19));
    EXPECT_EQ(dw.componentOfNode(3), dw.componentOfNode(7));

    g.addEdge(13, 15, 0);
    dw.update(GraphEvent(GraphEvent::EDGE_ADDITION, 13, 15, 0));
    EXPECT_EQ(4, dw.numberOfComponents());
    EXPECT_EQ(dw.componentOfNode(14), dw.componentOfNode(15));
    EXPECT_NE(dw.componentOfNode(15), dw.componentOfNode(0));

    // Create batch and update
    std::vector<GraphEvent> batch {
        GraphEvent(GraphEvent::EDGE_ADDITION, 15, 19, 0),
        GraphEvent(GraphEvent::EDGE_REMOVAL, 6, 7, 0),
        GraphEvent(GraphEvent::EDGE_ADDITION, 7, 17, 0),
        GraphEvent(GraphEvent::EDGE_REMOVAL, 2, 4, 0)
    };

    for (auto e : batch) {
        if (e.type == GraphEvent::EDGE_ADDITION){
            g.addEdge(e.u, e.v);
        }
        if (e.type == GraphEvent::EDGE_REMOVAL) {
            g.removeEdge(e.u, e.v);
        }
    }

    dw.updateBatch(batch);
    EXPECT_EQ(4, dw.numberOfComponents());
    EXPECT_TRUE(dw.componentOfNode(0) != dw.componentOfNode(14));
    EXPECT_TRUE(dw.componentOfNode(9) == dw.componentOfNode(11));
    EXPECT_TRUE(dw.componentOfNode(3) != dw.componentOfNode(9));
}


TEST_F(ConnectedComponentsGTest, testDynWeaklyConnectedComponents) {
    // Read graph
    KONECTGraphReader reader;
    Graph G = reader.read("input/foodweb-baydry.konect");
    DynWeaklyConnectedComponents dw(G);
    dw.run();
    WeaklyConnectedComponents wc(G);
    wc.run();

    DEBUG("Number of components: ", dw.numberOfComponents());
    // Testing static run.
    EXPECT_EQ(wc.numberOfComponents(), dw.numberOfComponents());

    int numberOfTests = 1000;
    // Probability to perform an edge insertion or removal.
    float p = 0.5;
    for (int i = 0; i < numberOfTests; ++i) {
        node u = 0;
        node v = 1;
        // Perform edge insertion
        if (((double) rand() / (RAND_MAX)) > p) {
            while (G.hasEdge(u, v)) {
                u = GraphTools::randomNode(G);
                v = GraphTools::randomNode(G);
            }
            G.addEdge(u, v);
            dw.update(GraphEvent(GraphEvent::EDGE_ADDITION, u, v, 0));
        }
        else {
            while (!G.hasEdge(u, v) || u == v) {
                std::pair<node, node> edge = GraphTools::randomEdge(G);
                u = edge.first;
                v = edge.second;
            }
            G.removeEdge(u, v);
            dw.update(GraphEvent(GraphEvent::EDGE_REMOVAL, u, v, 0));
        }
    }
    wc.run();
    EXPECT_EQ(wc.numberOfComponents(), dw.numberOfComponents());

    // Testing batch update.
    std::vector<GraphEvent> batch(numberOfTests);
    for (int i = 0; i < numberOfTests; ++i) {
        node u = GraphTools::randomNode(G);
        node v = GraphTools::randomNode(G);
        if (((double) rand() / (RAND_MAX)) > p) {
            while (G.hasEdge(u, v) || u == v) {
                u = GraphTools::randomNode(G);
                v = GraphTools::randomNode(G);
            }
            batch[i] = GraphEvent(GraphEvent::EDGE_ADDITION, u, v, 0);
            G.addEdge(u, v);
        }
        else {
            while (!G.hasEdge(u, v) || u == v) {
                std::pair<node, node> edge = GraphTools::randomEdge(G);
                u = edge.first;
                v = edge.second;
            }
            G.removeEdge(u, v);
            batch[i] = GraphEvent(GraphEvent::EDGE_REMOVAL, u, v, 0);
        }
    }

    dw.updateBatch(batch);
    wc.run();
    EXPECT_EQ(wc.numberOfComponents(), dw.numberOfComponents());
}

TEST_F(ConnectedComponentsGTest, testExtractLargestConnectedComponent) {
    Graph G(8);

    G.addEdge(0, 1);
    G.addEdge(2, 1);
    G.addEdge(3, 1);
    G.addEdge(4, 1);

    G.addEdge(5, 6);
    Graph G1(G);

    auto lcc = ConnectedComponents::extractLargestConnectedComponent(G, true);
    EXPECT_EQ(lcc.numberOfNodes(), 5);
    EXPECT_EQ(lcc.upperNodeIdBound(), 5);
    EXPECT_EQ(lcc.numberOfEdges(), 4);

    G.removeNode(0);
    lcc = ConnectedComponents::extractLargestConnectedComponent(G, false);
    lcc = ConnectedComponents::extractLargestConnectedComponent(lcc, true);
    EXPECT_EQ(lcc.numberOfNodes(), 4);
    EXPECT_EQ(lcc.upperNodeIdBound(), 4);
    EXPECT_EQ(lcc.numberOfEdges(), 3);
    node u = 0;
    lcc.forNodes([&u](const node v) { EXPECT_EQ(u++, v); });

    G1 = ConnectedComponents::extractLargestConnectedComponent(G1, false);
    EXPECT_EQ(G1.numberOfNodes(), 5);
    EXPECT_EQ(G1.upperNodeIdBound(), 8);
    EXPECT_EQ(G1.numberOfEdges(), 4);
}

TEST_F(ConnectedComponentsGTest, testExtractLargestStronglyConnectedComponent) {
    Aux::Random::setSeed(42, true);
    for (double p : {0.02, 0.025, 0.03}) {
        const auto G = ErdosRenyiGenerator(100, p, true).generate();
        auto lSCC = StronglyConnectedComponents::extractLargestStronglyConnectedComponent(G);
        StronglyConnectedComponents scc(G);
        scc.run();
        const auto sizes = scc.getPartition().subsetSizes();
        const count maxSize = *std::max_element(sizes.begin(), sizes.end());

        EXPECT_EQ(lSCC.numberOfNodes(), maxSize);
        StronglyConnectedComponents sccCheck(lSCC);
        sccCheck.run();
        EXPECT_EQ(sccCheck.numberOfComponents(), 1);

        lSCC = StronglyConnectedComponents::extractLargestStronglyConnectedComponent(G, true);
        node u = 0;
        lSCC.forNodes([&u](node v) {
            EXPECT_EQ(u, v);
            ++u;
        });
    }
}

} /* namespace NetworKit */
