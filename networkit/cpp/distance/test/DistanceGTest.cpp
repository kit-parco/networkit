// no-networkit-format
/*
 * DistanceGTest.cpp
 *
 *  Created on: Sep 04, 2015
 *      Author: Maximilian Vogel
 */

#include <gtest/gtest.h>

#include <networkit/distance/APSP.hpp>
#include <networkit/distance/AStar.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/BidirectionalBFS.hpp>
#include <networkit/distance/BidirectionalDijkstra.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/EffectiveDiameter.hpp>
#include <networkit/distance/EffectiveDiameterApproximation.hpp>
#include <networkit/distance/HopPlotApproximation.hpp>
#include <networkit/distance/IncompleteDijkstra.hpp>
#include <networkit/distance/NeighborhoodFunction.hpp>
#include <networkit/distance/NeighborhoodFunctionApproximation.hpp>
#include <networkit/distance/NeighborhoodFunctionHeuristic.hpp>
#include <networkit/distance/SPSP.hpp>

#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {
class DistanceGTest: public testing::Test {};

TEST_F(DistanceGTest, testVertexDiameterPedantically) {
    DorogovtsevMendesGenerator generator(1000);
    Graph G1 = generator.generate();
    Graph G = Graph(G1, true, false);
    Diameter diam(G, DiameterAlgo::estimatedPedantic);
    diam.run();
    count vd = diam.getDiameter().first;
    EXPECT_EQ(1000, vd);
}

TEST_F(DistanceGTest, testAStar) {

    // Builds a mesh graph with the given number of rows and columns
    auto buildMesh = [](count rows, count cols) -> Graph {
        Graph G(rows * cols, false, false);

        for (count i = 0; i < rows; ++i) {
            for (count j = 0; j < cols; ++j) {
                if (j < cols - 1) {
                    G.addEdge(i * cols + j, i * cols + j + 1);
                }
                if (i < rows - 1) {
                    G.addEdge(i * cols + j, (i + 1) * cols + j);
                }
            }
        }

        return G;
    };

    // Test AStar on a mesh for every pair of nodes (u, v), u != v.
    auto testMesh = [&](count rows, count cols) {
        auto G = buildMesh(rows, cols);

        // Test AStar for every pair.
        G.forNodePairs([&](node source, node target) {

            // Three distance heuristics
            std::vector<double> zeroDist(rows*cols, 0);
            std::vector<double> exactDist(rows*cols);
            std::vector<double> eucledianDist(rows*cols);
            for (node u = 0; u < rows*cols; ++u) {
                count rowU = u / cols;
                count colU = u % cols;
                count rowT = target / cols;
                count colT = target % cols;
                exactDist[u] = static_cast<double>((rowU > rowT ? rowU - rowT : rowT - rowU) +
                       (colU > colT ? colU - colT : colT - colU));
                double rowDiff = std::abs(static_cast<double>(u / cols) - static_cast<double>(rowT));
                double colDiff = std::abs(static_cast<double>(u % cols) - static_cast<double>(rowT));
                eucledianDist[u] = sqrt(std::pow(rowDiff, 2) + std::pow(colDiff, 2));
            };
            BFS bfs(G, source, true, false, target);
            bfs.run();

            std::vector<std::vector<double>> heuristics;

            for (auto heu : heuristics) {
                AStar astar(G, heu, source, target, true);
                astar.run();
                EXPECT_EQ(astar.getDistance(), bfs.distance(target));

                auto path = astar.getPath();
                EXPECT_EQ(path.size(), bfs.getPath(target).size() - 2);

                if (!path.size()) {
                    return;
                }
                for (size_t i = 0; i < path.size() - 1; ++i) {
                    EXPECT_TRUE(G.hasEdge(path[i], path[i + 1]));
                }
            }
        });
    };

    // Test cases
    testMesh(10, 10);
    testMesh(1, 15);
    testMesh(25, 5);
}

TEST_F(DistanceGTest, testIncompleteDijkstra) {
    Aux::Random::setSeed(42, false);
    for (auto directed : {true, false}) {
        for (auto weighted : {true, false}) {
            auto G = ErdosRenyiGenerator(500, 0.05, directed).generate();
            if (weighted) {
                G = GraphTools::toWeighted(G);
                G.forEdges([&G](const node u, const node v) {
                    G.setWeight(u, v, Aux::Random::probability());
                });
            }

            G.forNodes([&](const node source) {
                Dijkstra dij(G, source, false, false);
                dij.run();
                const auto dists = dij.getDistances();
                const count reachable = std::count_if(dists.begin(), dists.end(),
                                                      [](const edgeweight dist) {
                    return dist != std::numeric_limits<edgeweight>::max();
                });

                const std::vector<node> sources({source});
                IncompleteDijkstra iDij(&G, sources);

                for (count i = 0; i < reachable; ++i) {
                    EXPECT_TRUE(iDij.hasNext());
                    const auto next = iDij.next();
                    EXPECT_DOUBLE_EQ(next.second, dists[next.first]);
                }

                EXPECT_FALSE(iDij.hasNext());
            });
        }
    }
}

TEST_F(DistanceGTest, testBidirectionalBFS) {
    Aux::Random::setSeed(42, false);
    Graph G = ErdosRenyiGenerator(500, 0.02, false).generate();
    Graph G1 = ErdosRenyiGenerator(500, 0.05, true).generate();
    auto testGraph = [&](const Graph &G) {
        node source = GraphTools::randomNode(G);
        node target = GraphTools::randomNode(G);
        while (source == target)
            target = GraphTools::randomNode(G);
        BFS bfs(G, source, true, false, target);
        bfs.run();
        BidirectionalBFS bbfs(G, source, target, true);
        bbfs.run();
        if (bfs.distance(target) < G.upperNodeIdBound())
            EXPECT_EQ(bfs.distance(target), bbfs.getDistance());
        else
            EXPECT_EQ(std::numeric_limits<count>::max(), bbfs.getDistance());
        auto path = bbfs.getPath();
        EXPECT_EQ(path.size(), bbfs.getDistance() - 1);
        if (path.size()) {
            EXPECT_TRUE(G.hasEdge(source, path.front()));
            EXPECT_TRUE(G.hasEdge(path.back(), target));
            for (count i = 1; i < path.size()-1; ++i)
                EXPECT_TRUE(G.hasEdge(path[i], path[i + 1]));
        } else
            EXPECT_TRUE(G.hasEdge(source, target));
    };
    for (int i = 0; i < 10; ++i) {
        testGraph(G);
        testGraph(G1);
    }
}

TEST_F(DistanceGTest, testBidirectionalDijkstra) {
    Aux::Random::setSeed(42, false);
    count n = 300;
    Graph G(n, true, false);
    Graph G1(n, true, true);
    Graph gen = ErdosRenyiGenerator(n, 0.05, false).generate();
    Graph gen1 = ErdosRenyiGenerator(n, 0.05, true).generate();
    gen.forEdges([&](node u, node v) {
        G.setWeight(u, v, Aux::Random::probability());
    });
    gen1.forEdges([&](node u, node v) {
        G1.setWeight(u, v, Aux::Random::probability());
    });

    auto testGraph = [&](const Graph &G) {
        node source = GraphTools::randomNode(G);
        node target = GraphTools::randomNode(G);
        while (source == target)
            target = GraphTools::randomNode(G);
        BidirectionalDijkstra bdij(G, source, target, true);
        bdij.run();
        Dijkstra dij(G, source, true, true, target);
        dij.run();

        EXPECT_NEAR(bdij.getDistance(), dij.distance(target), 1e-6);
        auto path = bdij.getPath();

        if (path.size()) {
            EXPECT_TRUE(G.hasEdge(source, path.front()));
            for (size_t i = 0; i < path.size() - 1; ++i)
                EXPECT_TRUE(G.hasEdge(path[i], path[i + 1]));
            EXPECT_TRUE(G.hasEdge(path.back(), target));
        } else
            EXPECT_TRUE(G.hasEdge(source, target));
    };

    for (int i = 0; i < 10; ++i) {
        testGraph(G);
        testGraph(G1);
    }
}

TEST_F(DistanceGTest, testExactDiameter) {
    using namespace std;

    vector<pair<string, count>> testInstances= {pair<string, count>("lesmis", 14),
                                               pair<string, count>("jazz", 6),
                                               pair<string, count>("celegans_metabolic", 7)
                                              };

    for (auto testInstance : testInstances) {
        METISGraphReader reader;
        Graph G = reader.read("input/" + testInstance.first + ".graph");
        Diameter diam(G, DiameterAlgo::exact);
        diam.run();
        count diameter = diam.getDiameter().first;
        EXPECT_EQ(diameter, testInstance.second);
    }
}


TEST_F(DistanceGTest, testEstimatedDiameterRange) {
    using namespace std;

   vector<pair<string, count>> testInstances= {
                                               pair<string, count>("celegans_metabolic", 7),
                                               pair<string, count>("jazz", 6)
                                              };

    for (auto testInstance : testInstances) {
        METISGraphReader reader;
        Graph G = reader.read("input/" + testInstance.first + ".graph");
        Diameter diam(G, DiameterAlgo::estimatedRange, 0.1);
        diam.run();
        std::pair<count, count> range = diam.getDiameter();
        EXPECT_GE(testInstance.second, range.first);
        EXPECT_LE(testInstance.second, range.second);
    }
}
TEST_F(DistanceGTest, testPedanticDiameterErdos) {
    count n = 5000;
    ErdosRenyiGenerator gen(n,0.001);
    Graph G1 = gen.generate();
    Diameter diam(G1, DiameterAlgo::estimatedPedantic);
    diam.run();
    count diameter = diam.getDiameter().first;
    ASSERT_LE(diameter, n);
}


TEST_F(DistanceGTest, testEffectiveDiameterMinimal) {
    // Minimal example from the paper
    Graph G(5);
    G.addEdge(0,1);
    G.addEdge(1,2);
    G.addEdge(2,3);
    G.addEdge(3,4);
    G.addEdge(4,0);
    EffectiveDiameterApproximation aef(G);
    aef.run();
    double effective = aef.getEffectiveDiameter();
    Diameter diam(G, DiameterAlgo::exact);
    diam.run();
    count exact = diam.getDiameter().first;
    EXPECT_LE(effective, exact);
}

TEST_F(DistanceGTest, testEffectiveDiameter) {

using namespace std;

vector<string> testInstances= {"celegans_metabolic", "jazz", "lesmis"};

for (auto testInstance : testInstances) {
    METISGraphReader reader;
    Graph G = reader.read("input/" + testInstance + ".graph");
    EffectiveDiameterApproximation aef(G);
    aef.run();
    double effective = aef.getEffectiveDiameter();
    Diameter diam(G, DiameterAlgo::exact);
    diam.run();
    count exact = diam.getDiameter().first;
    EXPECT_LE(effective, exact);
}
}

TEST_F(DistanceGTest, testEffectiveDiameterExact) {

    using namespace std;

    vector<string> testInstances= {"celegans_metabolic", "jazz", "lesmis"};

    for (auto testInstance : testInstances) {
        METISGraphReader reader;
        Graph G = reader.read("input/" + testInstance + ".graph");
        EffectiveDiameter ed(G);
        ed.run();
        double effective = ed.getEffectiveDiameter();
        Diameter diam(G, DiameterAlgo::exact);
        diam.run();
        count exact = diam.getDiameter().first;
        EXPECT_LE(effective, exact);
    }

    const double tol = 1e-3;

    /* Graph: n=20, threshold: 20*0.9 = 18 nodes
        1--3--5--7---9
        |  |  |  |   |
        2--4--6--8--10
            |     |
            11----12
                |
            13--14--15
                |
            18--16--17--19
                    |
                    20
    Number of steps needed per node: (1-20)
    (7+6+6+5+6+5+5+4+6+5+4+4+5+4+5+5+6+6+7+7) / 20 = 5.4
    */
        count n1 = 20;
        Graph G1(n1);

        G1.addEdge(0,1);
        G1.addEdge(0,2);
        G1.addEdge(1,3);
        G1.addEdge(2,3);
        G1.addEdge(2,4);
        G1.addEdge(3,5);
        G1.addEdge(3,10);
        G1.addEdge(4,5);
        G1.addEdge(4,6);
        G1.addEdge(5,7);
        G1.addEdge(6,8);
        G1.addEdge(6,7);
        G1.addEdge(7,9);
        G1.addEdge(7,11);
        G1.addEdge(8,9);
        G1.addEdge(10,11);
        G1.addEdge(11,13);
        G1.addEdge(12,13);
        G1.addEdge(13,14);
        G1.addEdge(13,15);
        G1.addEdge(15,16);
        G1.addEdge(15,17);
        G1.addEdge(16,18);
        G1.addEdge(16,19);

        EffectiveDiameter ed(G1);
        ed.run();
        double effective1 = ed.getEffectiveDiameter();
        EXPECT_NEAR(5.4, effective1, tol);

        /* Graph: n=21, threshold: 21*0.9 = 18.9 => 19 nodes
                    13---------------3
                        |               |
                    ---14--12--|        |
                    |   |   |  |        |
        1--21--18--16--15   |  |        |
                    |       |  |        |
            20--17------10--8        |
                    |       |  |        |
                19       9--7--5--6--4--11
                                        |
                                        2
    Number of steps needed per node: (1-21)
    (8+7+5+6+6+6+5+5+5+5+7+5+4+4+5+5+5+6+6+6+7) / 21 = 5.619047
    */
        count n2 = 21;
        Graph G2(n2);

        G2.addEdge(0,20);
        G2.addEdge(1,3);
        G2.addEdge(2,3);
        G2.addEdge(2,12);
        G2.addEdge(3,5);
        G2.addEdge(3,10);
        G2.addEdge(4,5);
        G2.addEdge(4,6);
        G2.addEdge(6,7);
        G2.addEdge(6,8);
        G2.addEdge(7,9);
        G2.addEdge(7,11);
        G2.addEdge(8,9);
        G2.addEdge(9,11);
        G2.addEdge(9,16);
        G2.addEdge(11,13);
        G2.addEdge(12,13);
        G2.addEdge(13,14);
        G2.addEdge(13,15);
        G2.addEdge(14,15);
        G2.addEdge(15,16);
        G2.addEdge(15,17);
        G2.addEdge(16,18);
        G2.addEdge(16,19);
        G2.addEdge(17,20);

        EffectiveDiameter ed2(G2);
        ed2.run();
        double effective2 = ed2.getEffectiveDiameter();
        EXPECT_NEAR(5.619047, effective2, tol);
}

TEST_F(DistanceGTest, testHopPlotApproximation) {
    using namespace std;

    vector<string> testInstances= {"celegans_metabolic", "lesmis"};

    const double tol = 1e-2;

    for (auto& testInstance : testInstances) {
        METISGraphReader reader;
        Graph G = reader.read("input/" + testInstance + ".graph");
        HopPlotApproximation hp(G);
        hp.run();
        map<count, double> hopPlot = hp.getHopPlot();
        for (count i=1; i < hopPlot.size(); i++) {
            EXPECT_LE(hopPlot[i-1], hopPlot[i]+tol);
        }
    }
}

TEST_F(DistanceGTest, testNeighborhoodFunctionApproximation) {
    METISGraphReader reader;
    Graph G = GraphTools::toUnweighted(reader.read("input/lesmis.graph"));
    NeighborhoodFunction nf(G);
    nf.run();
    auto exact = nf.getNeighborhoodFunction();
    NeighborhoodFunctionApproximation anf(G);
    anf.run();
    auto approximated = anf.getNeighborhoodFunction();
    EXPECT_EQ(exact.size(), approximated.size());
}

TEST_F(DistanceGTest, testNeighborhoodFunctionHeuristic) {
    METISGraphReader reader;
    Graph G = GraphTools::toUnweighted(reader.read("input/lesmis.graph"));
    NeighborhoodFunction nf(G);
    nf.run();
    auto exact = nf.getNeighborhoodFunction();
    NeighborhoodFunctionHeuristic anf(G);
    anf.run();
    auto heuristic = anf.getNeighborhoodFunction();
    EXPECT_EQ(exact.size(), heuristic.size());
}

TEST_F(DistanceGTest, testSPSP) {
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);
        for (bool directed : {false, true}) {
            for (bool weighted : {false, true}) {
                auto G = ErdosRenyiGenerator(100, 0.15, directed).generate();
                if (weighted) {
                    G = GraphTools::toWeighted(G);
                    G.forEdges([&](node u, node v) {
                        G.setWeight(u, v, Aux::Random::probability());
                    });
                }

                APSP apsp(G);
                apsp.run();
                const auto gt = apsp.getDistances();

                for (count nSources : {1, 10, 50, 100}) {
                    std::unordered_set<node> sources;
                    do {
                        sources.insert(GraphTools::randomNode(G));
                    } while (sources.size() < nSources);

                    SPSP spsp(G, sources.begin(), sources.end());
                    spsp.run();
                    const auto nodemap = spsp.getSourceIndexMap();
                    const auto distances = spsp.getDistances();

                    for (node source : sources) {
                        const auto &curDist = distances[nodemap.at(source)];
                        G.forNodes([&](node target) {
                            EXPECT_DOUBLE_EQ(gt[source][target], spsp.getDistance(source, target));
                            EXPECT_DOUBLE_EQ(gt[source][target], curDist[target]);
                        });
                    }
                }
            }
        }
    }
}

} /* namespace NetworKit */
