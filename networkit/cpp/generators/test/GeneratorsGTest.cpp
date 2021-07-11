// no-networkit-format
/*
 *  GeneratorsTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include <gtest/gtest.h>

#include <numeric>
#include <cmath>

#include <networkit/generators/ClusteredRandomGraphGenerator.hpp>
#include <networkit/generators/DynamicGraphSource.hpp>
#include <networkit/generators/DynamicBarabasiAlbertGenerator.hpp>
#include <networkit/generators/PubWebGenerator.hpp>
#include <networkit/generators/DynamicPubWebGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/generators/ChungLuGenerator.hpp>
#include <networkit/generators/HavelHakimiGenerator.hpp>
#include <networkit/generators/RmatGenerator.hpp>
#include <networkit/generators/BarabasiAlbertGenerator.hpp>
#include <networkit/generators/DynamicPathGenerator.hpp>
#include <networkit/generators/DynamicForestFireGenerator.hpp>
#include <networkit/generators/DynamicDorogovtsevMendesGenerator.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/WattsStrogatzGenerator.hpp>
#include <networkit/generators/RegularRingLatticeGenerator.hpp>
#include <networkit/generators/StochasticBlockmodel.hpp>
#include <networkit/generators/EdgeSwitchingMarkovChainGenerator.hpp>
#include <networkit/generators/LFRGenerator.hpp>
#include <networkit/generators/MocnikGenerator.hpp>
#include <networkit/generators/MocnikGeneratorBasic.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/generators/DynamicHyperbolicGenerator.hpp>

#include <networkit/viz/PostscriptWriter.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/community/PLP.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/io/METISGraphWriter.hpp>
#include <networkit/io/DotGraphWriter.hpp>
#include <networkit/io/GraphIO.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/dynamics/GraphUpdater.hpp>
#include <networkit/auxiliary/MissingMath.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/global/ClusteringCoefficient.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/community/Modularity.hpp>


namespace NetworKit {

class GeneratorsGTest: public testing::Test {
public:
    vector<double> getAngles(const DynamicHyperbolicGenerator &dynGen) { return dynGen.angles; }

    vector<double> getRadii(const DynamicHyperbolicGenerator &dynGen) { return dynGen.radii; }
};

TEST_F(GeneratorsGTest, testClusteredRandomGraphGenerator) {
    Aux::Random::setSeed(42, false);
    const count n = 100, c = 10;
    const double pin = 0.5, pout = 0.01;
    ClusteredRandomGraphGenerator gen(n, c, pin, pout);
    Graph G = gen.generate();
    Partition part = gen.getCommunities();
    count nCommunities = part.getSubsetIds().size();
    EXPECT_EQ(n, G.numberOfNodes());
    EXPECT_EQ(n, G.upperNodeIdBound());
    EXPECT_GE(nCommunities, 1);
    EXPECT_LE(nCommunities, c);
}

TEST_F(GeneratorsGTest, testClusteredRandomGraphGeneratorCompleteIsolatedCommunities) {
    const count n = 100, c = 10;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        ClusteredRandomGraphGenerator gen(n, c, 1, 0);
        auto G = gen.generate();
        EXPECT_EQ(n, G.numberOfNodes());
        EXPECT_EQ(n, G.upperNodeIdBound());
        EXPECT_EQ(G.numberOfSelfLoops(), 0);
        {   // Check no multiple edges
            auto G1 = G;
            G1.removeMultiEdges();
            EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        }

        const auto commPartition = gen.getCommunities();
        const auto communities = commPartition.getSubsets();
        EXPECT_GE(communities.size(), 1);
        EXPECT_LE(communities.size(), c);

        // Check that each community is an isolated clique
        const auto subsetSizes = commPartition.subsetSizeMap();
        for (const auto &community : communities) {
            for (node u : community) {
                const index uIdx = commPartition.subsetOf(u);
                EXPECT_EQ(G.degree(u), subsetSizes.at(uIdx) - 1);
                G.forNeighborsOf(u, [&](node v) {
                    EXPECT_EQ(uIdx, commPartition.subsetOf(v));
                });
            }
        }
    }
}

TEST_F(GeneratorsGTest, testClusteredRandomGraphGeneratorCompleteCommunities) {
    const count n = 100, c = 10;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        ClusteredRandomGraphGenerator gen(n, c, 1, 0.01);
        auto G = gen.generate();
        EXPECT_EQ(n, G.numberOfNodes());
        EXPECT_EQ(n, G.upperNodeIdBound());
        EXPECT_EQ(G.numberOfSelfLoops(), 0);
        {   // Check no multiple edges
            auto G1 = G;
            G1.removeMultiEdges();
            EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        }

        const auto commPartition = gen.getCommunities();
        const auto communities = commPartition.getSubsets();
        EXPECT_GE(communities.size(), 1);
        EXPECT_LE(communities.size(), c);

        // Check that each community is a clique
        const auto subsetSizes = commPartition.subsetSizeMap();
        for (const auto &community : communities) {
            for (node u : community) {
                const index uIdx = commPartition.subsetOf(u);
                EXPECT_GE(G.degree(u), subsetSizes.at(uIdx) - 1);
                count leftToVisit = community.size() - 1;
                G.forNeighborsOf(u, [&](node v) {
                    if (uIdx == commPartition.subsetOf(v))
                        --leftToVisit;
                });
                EXPECT_EQ(leftToVisit, 0);
            }
        }
    }
}

TEST_F(GeneratorsGTest, testClusteredRandomGraphGeneratorIndependentCommunities) {
    const count n = 100, c = 10;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        ClusteredRandomGraphGenerator gen(n, c, 0, 1);
        auto G = gen.generate();
        EXPECT_EQ(n, G.numberOfNodes());
        EXPECT_EQ(n, G.upperNodeIdBound());
        EXPECT_EQ(G.numberOfSelfLoops(), 0);
        {   // Check no multiple edges
            auto G1 = G;
            G1.removeMultiEdges();
            EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        }

        const auto commPartition = gen.getCommunities();
        const auto communities = commPartition.getSubsets();
        EXPECT_GE(communities.size(), 1);
        EXPECT_LE(communities.size(), c);

        // Check that each community is an independent set, and the vertices in
        // each community are connected to all the vertices outside the community.
        const auto subsetSizes = commPartition.subsetSizeMap();
        for (const auto &community : communities) {
            for (node u : community) {
                EXPECT_EQ(G.degree(u), n - community.size());
                G.forNeighborsOf(u, [&](node v) {
                    EXPECT_NE(commPartition.subsetOf(u), commPartition.subsetOf(v));
                });
            }
        }
    }
}

TEST_F(GeneratorsGTest, testClusteredRandomGraphGeneratorSparseCommunities) {
    const count n = 100, c = 10;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        ClusteredRandomGraphGenerator gen(n, c, 0.01, 1);
        auto G = gen.generate();
        EXPECT_EQ(n, G.numberOfNodes());
        EXPECT_EQ(n, G.upperNodeIdBound());
        EXPECT_EQ(G.numberOfSelfLoops(), 0);
        {   // Check no multiple edges
            auto G1 = G;
            G1.removeMultiEdges();
            EXPECT_EQ(G.numberOfEdges(), G1.numberOfEdges());
        }

        const auto commPartition = gen.getCommunities();
        const auto communities = commPartition.getSubsets();
        EXPECT_GE(communities.size(), 1);
        EXPECT_LE(communities.size(), c);

        // Check that the vertices in each community are connected to all the
        // vertices outside the community.
        const auto subsetSizes = commPartition.subsetSizeMap();
        for (const auto &community : communities) {
            for (node u : community) {
                EXPECT_GE(G.degree(u), n - community.size());
                count leftToVisit = n - community.size();
                G.forNeighborsOf(u, [&](node v) {
                    if (commPartition.subsetOf(u) != commPartition.subsetOf(v))
                        --leftToVisit;
                });
                EXPECT_EQ(leftToVisit, 0);
            }
        }
    }
}


TEST_F(GeneratorsGTest, testDynamicBarabasiAlbertGeneratorSingleStep) {
    count k = 2; // number of edges added per node
    DynamicGraphSource* gen = new DynamicBarabasiAlbertGenerator(k);
    GraphEventProxy* Gproxy = gen->newGraph();
    Graph* G = Gproxy->G;

    gen->initializeGraph();

    count nPre = G->numberOfNodes();
    count mPre = G->numberOfEdges();
    EXPECT_EQ(k, nPre) << "graph should have been initialized to k nodes";
    EXPECT_EQ(k - 1, mPre) << "graph should have been initialized to a path of k nodes which means k-1 edges";

    // perform single preferential attachment step
    gen->generate();

    count nPost = G->numberOfNodes();
    count mPost = G->numberOfEdges();
    EXPECT_EQ(nPre + 1, nPost) << "one more node should have been added";
    EXPECT_EQ(mPre + k, mPost) << "k edges should have been added";

    delete gen;
    delete Gproxy;
    delete G;
}

TEST_F(GeneratorsGTest, testDynamicBarabasiAlbertGenerator) {
    DynamicGraphSource* gen = new DynamicBarabasiAlbertGenerator(2);

    GraphEventProxy* Gproxy = gen->newGraph();
    Graph* G = Gproxy->G;

    gen->initializeGraph();

    EXPECT_EQ(2u, G->numberOfNodes()) << "initially the generator creates two connected nodes";
    EXPECT_EQ(1u, G->numberOfEdges()) << "initially the generator creates two connected nodes";

    count n = 100;

    gen->generateWhile([&]() {
                return ( G->numberOfNodes() < n );
            });

    EXPECT_EQ(n, G->numberOfNodes());
    DEBUG("m = " , G->numberOfEdges());

    // resume generator

    gen->generateWhile([&]() {
        return (G->numberOfNodes() < 2 * n);
    });
    EXPECT_EQ(2 * n, G->numberOfNodes());

    delete gen;
    delete Gproxy;
    delete G;
}


TEST_F(GeneratorsGTest, viewDynamicBarabasiAlbertGenerator) {
    DynamicGraphSource* gen = new DynamicBarabasiAlbertGenerator(2);
    GraphEventProxy* Gproxy = gen->newGraph();
    Graph* G = Gproxy->G;
    gen->initializeGraph();
    count n = 42;
    gen->generateWhile([&]() {
                return ( G->numberOfNodes() < n );
            });
    METISGraphWriter writer;
    writer.write(*G, "output/BATest.graph");

    delete gen;
    delete Gproxy;
    delete G;
}

TEST_F(GeneratorsGTest, testStaticPubWebGenerator) {
    Aux::Random::setSeed(42, false);

    count n = 450;
    count numCluster = 9;
    count maxNumNeighbors = 36;
    float rad = 0.075;

    PubWebGenerator gen(n, numCluster, rad, maxNumNeighbors);
    Graph G = gen.generate();
    auto coordinates = gen.moveCoordinates();
    EXPECT_EQ(n, G.numberOfNodes()) << "number of generated nodes";

    // check degree
    G.forNodes([&](node v) {
        EXPECT_LE(G.degree(v), maxNumNeighbors) << "maximum degree";
    });

    // 1-clustering
    ClusteringGenerator clusterGen;
    Partition oneClustering = clusterGen.makeOneClustering(G);
    EXPECT_EQ(G.numberOfNodes(),oneClustering.numberOfElements());

    // output to EPS file
    PostscriptWriter psWriter(true);
    psWriter.write(G, coordinates, oneClustering, "output/pubweb.eps");

    // clustering
    PLM clusterAlgo(G, false, 1.0, "none randomized");
    clusterAlgo.run();
    Partition clustering = clusterAlgo.getPartition();
    EXPECT_EQ(G.numberOfNodes(),clustering.numberOfElements());
    psWriter.write(G, coordinates, clustering, "output/pubweb-clustered-PLM.eps");

    Modularity mod;
    double modVal = mod.getQuality(clustering, G);
    EXPECT_GE(modVal, 0.2) << "modularity of clustering";
    DEBUG("Modularity of clustering: " , modVal);
    DEBUG("Total edge weight: " , G.totalEdgeWeight());
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testDynamicPubWebGenerator) {
    count nSteps = 5;
    count n = 200;
    count numCluster = 10;
    count maxNumNeighbors = 40;
    float rad = 0.08;

    DynamicPubWebGenerator dynGen(n, numCluster, rad, maxNumNeighbors, false);
    Graph G = dynGen.getGraph();
    GraphUpdater gu(G);
    std::vector<GraphEvent> stream;

    for (index i = 1; i <= nSteps; ++i) {
        stream = dynGen.generate(1);
        DEBUG("updating graph");
        gu.update(stream);

        DEBUG("updated graph, new (n, m) = (" , G.numberOfNodes() , ", " , G.numberOfEdges() , ")");
        edgeweight tew = G.totalEdgeWeight();
        DEBUG("1/2 graph volume: ", tew);
        EXPECT_GT(tew, 0);

        // update coordinates
        #ifndef NETWORKIT_RELEASE_LOGGING
        {
            PostscriptWriter psWriter(true);
            std::stringstream ss;
            ss << "output/pubweb-" << std::setw(4) << std::setfill('0') << i <<  ".eps";
            psWriter.write(G, dynGen.getCoordinates(), ss.str());
        }
        #endif
    }
}

/**
 * Testing the dynamic hyperbolic generator with fixed parameters and changing node positions
 */
TEST_F(GeneratorsGTest, testDynamicHyperbolicGeneratorOnMovedNodes) {
    Aux::Random::setSeed(0, false);

    //set up dynamic parameters
    int nSteps = 20;
    const count n = 500;
    const double k = 6;
    const double alpha = 1;
    //const double exp = 2*alpha+1;
    const double T = 0;
    const double R = HyperbolicSpace::getTargetRadius(n, n*k/2, alpha, T);

    double movedShare = 1;
    double moveDistance = 0.1;

    //set up initial node positions
    vector<double> angles(n, -1);
    vector<double> radii(n, -1);
    HyperbolicSpace::fillPoints(angles, radii, R, alpha);
    DynamicHyperbolicGenerator dynGen(angles, radii, R, alpha, T, movedShare, moveDistance);

    //generate starting graph
    Graph G = HyperbolicGenerator().generate(angles, radii, R);
    count initialEdgeCount = G.numberOfEdges();
    count expected = n*HyperbolicSpace::getExpectedDegree(n, alpha, R)*0.5;
    EXPECT_NEAR(initialEdgeCount, expected, expected/5);
    GraphUpdater gu(G);
    std::vector<GraphEvent> stream;

    for (int i = 0; i < nSteps; i++) {
        //move nodes and generate stream of affected edges
        stream = dynGen.generate(1);
        DEBUG("Edges: ", G.numberOfEdges());
        for (auto event : stream) {
            EXPECT_TRUE(event.type == GraphEvent::EDGE_REMOVAL || event.type == GraphEvent::EDGE_ADDITION || event.type == GraphEvent::TIME_STEP);
            if (event.type == GraphEvent::EDGE_REMOVAL) {
                EXPECT_TRUE(G.hasEdge(event.u, event.v));
            }
            //only present nodes can be affected, no new nodes are introduced
            if (event.type != GraphEvent::TIME_STEP){
                EXPECT_LT(event.u, G.upperNodeIdBound());
            }
        }
        gu.update(stream);
        EXPECT_TRUE(G.checkConsistency());
    }

    //update moved nodes
    angles = getAngles(dynGen);
    radii = getRadii(dynGen);
    Graph comparison = HyperbolicGenerator().generate(angles, radii, R);
    EXPECT_EQ(G.numberOfEdges(), comparison.numberOfEdges());

    //heuristic criterion: Number of edges may change, but should not change much
    EXPECT_NEAR(G.numberOfEdges(), initialEdgeCount, initialEdgeCount/5);
}

/**
 * creates a series of pictures visualizing the effect of the dynamic hyperbolic generator
 */
TEST_F(GeneratorsGTest, testDynamicHyperbolicVisualization) {
    count n = 300;
    count nSteps = 20;

    const double k = 6;
    const double alpha = 1;
    //const double exp = 2*alpha+1;
    const double T = 0;
    const double R = HyperbolicSpace::getTargetRadius(n, n*k/2, alpha, T);

    double movedShare = 0.2;
    double moveDistance = 1;
    vector<double> angles(n);
    vector<double> radii(n);

    HyperbolicSpace::fillPoints(angles, radii, R, alpha);

    DynamicHyperbolicGenerator dynGen(angles, radii, R, alpha, T, movedShare, moveDistance);
    Graph G = dynGen.getGraph();

    GraphUpdater gu(G);
    std::vector<GraphEvent> stream;

    for (index i = 0; i < nSteps; i++) {
        stream = dynGen.generate(1);
        DEBUG("Edges: ", G.numberOfEdges());
        for (auto event : stream) {
            EXPECT_TRUE(event.type == GraphEvent::EDGE_REMOVAL || event.type == GraphEvent::EDGE_ADDITION || event.type == GraphEvent::TIME_STEP);
        }
        gu.update(stream);

        auto coordinates = dynGen.getCoordinates();

        #ifndef NETWORKIT_RELEASE_LOGGING
        {
            std::stringstream ss;
            ss << "output/hyperbolic-" << std::setw(4) << std::setfill('0') << i << ".eps";

            PostscriptWriter psWriter(true);
            psWriter.write(G, coordinates, ss.str());
        }
        #endif
    }
}

TEST_F(GeneratorsGTest, testBarabasiAlbertGeneratorOriginal) {
    count k = 3;
    count nMax = 100;
    count n0 = 3;

    BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0, false);
    Graph G = BarabasiAlbert.generate();
    EXPECT_FALSE(G.isEmpty());

    EXPECT_EQ(nMax, G.numberOfNodes());
    EXPECT_EQ( ((n0-1) + ((nMax - n0) * k)), G.numberOfEdges());
    EXPECT_TRUE(G.checkConsistency());

    Graph initGraph(4);
    initGraph.addEdge(0,1);
    initGraph.addEdge(2,1);
    initGraph.addEdge(2,3);
    initGraph.addEdge(0,3);
    BarabasiAlbert = BarabasiAlbertGenerator(k, nMax, initGraph, false);
    G = BarabasiAlbert.generate();

    EXPECT_EQ(nMax, G.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), (nMax - initGraph.numberOfNodes()) * k + initGraph.numberOfEdges());
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testBarabasiAlbertGeneratorConstructor) {
    // k > nMax
    EXPECT_THROW(BarabasiAlbertGenerator generator(10, 9, 8, false), std::runtime_error);
    EXPECT_THROW(BarabasiAlbertGenerator generator(10, 9, 8, true), std::runtime_error);

    // n0 > nMax
    EXPECT_THROW(BarabasiAlbertGenerator generator(5, 9, 10, false), std::runtime_error);
    EXPECT_THROW(BarabasiAlbertGenerator generator(5, 9, 10, true), std::runtime_error);

    // n0 = initGraph.numberOfNodes() > nMax
    Graph initGraph(10);
    EXPECT_THROW(BarabasiAlbertGenerator generator(6, 9, initGraph, false), std::runtime_error);
    EXPECT_THROW(BarabasiAlbertGenerator generator(6, 9, initGraph, true), std::runtime_error);

    // initGraph, k > nMax
    initGraph = Graph(6);
    EXPECT_THROW(BarabasiAlbertGenerator generator(10, 9, initGraph, false), std::runtime_error);
    EXPECT_THROW(BarabasiAlbertGenerator generator(10, 9, initGraph, true), std::runtime_error);

    // initGraph, original method, initGraph.numberOfNodes() < k
    EXPECT_THROW(BarabasiAlbertGenerator generator(8, 9, initGraph, false), std::runtime_error);

    // initGraph does not have consecutive node ids
    initGraph.removeNode(0);
    EXPECT_THROW(BarabasiAlbertGenerator generator(3, 9, initGraph, false), std::runtime_error);
    EXPECT_THROW(BarabasiAlbertGenerator generator(3, 9, initGraph, false), std::runtime_error);
}

TEST_F(GeneratorsGTest, testBarabasiAlbertGeneratorBatagelj) {
    count k = 3;
    count nMax = 1000;
    count n0 = 3;

    BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0, true);
    Graph G = BarabasiAlbert.generate();

    EXPECT_EQ(nMax, G.numberOfNodes());
    EXPECT_LE(G.numberOfEdges(), nMax * k);
    EXPECT_TRUE(G.checkConsistency());

    Graph initGraph(4);
    initGraph.addEdge(0,1);
    initGraph.addEdge(2,1);
    initGraph.addEdge(2,3);
    initGraph.addEdge(0,3);
    BarabasiAlbert = BarabasiAlbertGenerator(k, nMax, initGraph, true);
    G = BarabasiAlbert.generate();

    EXPECT_EQ(nMax, G.numberOfNodes());
    EXPECT_LE(G.numberOfEdges(), nMax * k);
    EXPECT_TRUE(G.checkConsistency());

}

TEST_F(GeneratorsGTest, generatetBarabasiAlbertGeneratorGraph) {
        count k = 3;
        count nMax = 1000;
        count n0 = 3;

        BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);

        Graph G = BarabasiAlbert.generate();
        GraphIO io;
        io.writeAdjacencyList(G, "output/"
                "BarabasiGraph.txt");
}

TEST_F(GeneratorsGTest, testDynamicPathGenerator) {
    count nSteps = 42;
    DynamicPathGenerator gen;
    auto stream = gen.generate(nSteps);
    EXPECT_EQ(stream.size(),nSteps * 3 + 1);
    #ifndef NETWORKIT_RELEASE_LOGGING
        for (auto ev : stream) {
            TRACE(ev.toString());
        }
    #endif
}

TEST_F(GeneratorsGTest, testErdosRenyiGenerator) {
    count n = 2000;
    double p = 1.5 * (log(n) / (double) n);

    ErdosRenyiGenerator generator(n, p);
    Graph G = generator.generate();
    EXPECT_EQ(n, G.numberOfNodes());
    EXPECT_FALSE(G.isEmpty());
    EXPECT_TRUE(G.checkConsistency());

    count nPairs = (n * (n-1)) / 2;
    count nEdges = G.numberOfEdges();
    EXPECT_GE(nEdges, 0.75 * p * nPairs);
    EXPECT_LE(nEdges, 1.25 * p * nPairs);

    DEBUG("Number of edges with probability " , p , " (actual/expected): " , nEdges , " / " , (nPairs * p));
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testRmatGeneratorException) {
    count scale = 9;
    count edgeFactor = 12;
    double a = 0.51;
    double b = 0.12;
    double c = 0.12;
    double d = 0.2;

    EXPECT_THROW(RmatGenerator rmat(scale, edgeFactor, a, b, c, d), std::runtime_error);
}

TEST_F(GeneratorsGTest, testRmatGenerator) {
    count scale = 9;
    count n = (1 << scale);
    count edgeFactor = 12;
    double a = 0.51;
    double b = 0.12;
    double c = 0.12;
    double d = 0.25;

    RmatGenerator rmat(scale, edgeFactor, a, b, c, d);
    Graph G = rmat.generate();

    EXPECT_EQ(G.numberOfNodes(), n);
    EXPECT_LE(G.numberOfEdges(), n * edgeFactor);

    ClusteringCoefficient cc;
    double ccex = cc.exactGlobal(G);
    EXPECT_LE(ccex, 0.4);

    PLM clusterer(G, true);
    clusterer.run();
    Partition zeta = clusterer.getPartition();
    Modularity mod;
    double modVal = mod.getQuality(zeta, G);
    INFO("Modularity of R-MAT graph clustering: ", modVal);
    EXPECT_GE(modVal, 0.0);
    EXPECT_TRUE(G.checkConsistency());
}


TEST_F(GeneratorsGTest, testChungLuGenerator) {
    count n = 400;
    count maxDegree = n / 8;
    std::vector<count> sequence(n);
    count expVolume = 0;
    count actualVolume = 0;

    // fill sequence with random values (this is not power-law, of course!)
    for (index i = 0; i < n; ++i) {
        sequence[i] = rand() % maxDegree;
        expVolume += sequence[i];
    }

    ChungLuGenerator gen(sequence);
    Graph G = gen.generate();
    EXPECT_TRUE(G.checkConsistency());

    EXPECT_EQ(n, G.numberOfNodes());
    G.forNodes([&](node v) {
        actualVolume += G.degree(v);
    });

    INFO("expected volume: ", expVolume, ", actual volume: ", actualVolume);
}

TEST_F(GeneratorsGTest, testChungLuGeneratorDegreeConsistency) {
    count n = 1000;
    std::vector<count> vec;
    count maxDegree = n / 8;
    /* Creates a random sequence of weights */
    for (index i = 0; i < n; i++){
        int grad = Aux::Random::integer(1, maxDegree);
        vec.push_back(grad);
    }
    ChungLuGenerator generator(vec);
    Graph G = generator.generate();
    /* We check to see if the actual degrees of our nodes vary too much from the expected ones.
    * However, we need to sort the expected degrees first, since the algorithm does this as well
    * and the nodes with the highest degrees are added first. */
    Aux::Parallel::sort(vec.begin(), vec.end(), [](count a, count b){ return a > b;});
    /* Check if node degree is more than 50% off from the expected degree of that node. */
    // TODO Should we be looking for something better than a 50% range here?
    G.parallelForNodes([&] (node v) {
        EXPECT_NEAR(G.degree(v), vec[v], (0.5 * maxDegree));
    });
}

TEST_F(GeneratorsGTest, testChungLuGeneratorVolumeConsistency) {
    count n = 1000;
    std::vector<count> vec;
    count maxDegree = n / 8;
    count expectedVolume = 0;
    /* Creates a random sequence of weights */
    for (index i = 0; i < n; i++){
        int grad = Aux::Random::integer(1, maxDegree);
        vec.push_back(grad);
        expectedVolume += grad;
    }
    ChungLuGenerator generator(vec);
    Graph G = generator.generate();
    /* Check if volume is more than 10% off from the expected volume. */
    //TODO Is a 20% offset here sufficient? */
    EXPECT_NEAR(G.numberOfEdges() * 2, expectedVolume, 0.2 * expectedVolume);
}

TEST_F(GeneratorsGTest, testHavelHakimiGeneratorOnRandomSequence) {
    count n = 400;
    std::vector<count> sequence(n);
    bool realizable = false;

    std::mt19937_64 prng(1);
    std::uniform_int_distribution<count> deg_distr(0, n / 10 - 1);
    do {
        // fill sequence with random values (this is not power-law, of course!)
        std::generate(sequence.begin(), sequence.end(), [&] {return deg_distr(prng);});

        // check if sequence is realizable
        HavelHakimiGenerator hhgen(sequence);
        realizable = hhgen.isRealizable();

        if (realizable) {
            Graph G = hhgen.generate();
            EXPECT_TRUE(G.checkConsistency());
            count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
            EXPECT_EQ(volume, 2 * G.numberOfEdges());
        }
    } while (! realizable);
}

TEST_F(GeneratorsGTest, testHavelHakimiGeneratorOnRealSequence) {
    METISGraphReader reader;
    std::vector<std::string> graphs = {"input/jazz.graph",
            "input/lesmis.graph"}; //, "input/PGPgiantcompo.graph", "input/coAuthorsDBLP.graph"};

    for (auto path : graphs) {
        Graph G = reader.read(path);
        count n = G.numberOfNodes();
        std::vector<count> sequence(n);
        G.forNodes([&](node u){
            sequence[u] = G.degree(u);

        });

        HavelHakimiGenerator hhgen(sequence);
        Graph G2 = hhgen.generate();
        EXPECT_TRUE(G.checkConsistency());

        count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
        EXPECT_EQ(volume, 2 * G2.numberOfEdges());

        if (volume < 50000) {
            std::vector<count> testSequence(n);
            G2.forNodes([&](node u){
                testSequence[u] = G2.degree(u);
            });

            for (index i = 0; i < n; ++i) {
                EXPECT_EQ(sequence[i], testSequence[i]);
            }
        }
    }
}

TEST_F(GeneratorsGTest, testHavelHakimiGeneratorOnUnrealizableSequence) {
    std::vector<count> seq = {20, 10, 2, 2, 2, 2, 2, 2, 2, 2, 2};

    HavelHakimiGenerator hhgen(seq);
    EXPECT_THROW(hhgen.generate(), std::runtime_error);

    hhgen = HavelHakimiGenerator(seq, true);
    Graph G = hhgen.generate();

    G.forNodes([&](node u) {
        EXPECT_EQ(std::min<count>(seq[u], 10), G.degree(u));
    });
}


TEST_F(GeneratorsGTest, testDynamicForestFireGenerator) {
    Graph G1(0);
    GraphUpdater gu1(G1);
    std::vector<GraphEvent> stream;
    DynamicForestFireGenerator ffg1(0.0, false);
    stream = ffg1.generate(10);
    gu1.update(stream);
    EXPECT_TRUE(G1.checkConsistency());
    EXPECT_EQ(10u, G1.numberOfNodes());
    G1.forNodes([&](node u) {
        count c = 0;
        G1.forNeighborsOf(u, [&](node v) {
            if (v < u) {
                c += 1;
            }
        });
        if (u == 0) {
            EXPECT_EQ(0u, c);
        } else {
            EXPECT_EQ(1u, c);
        }
    });

    Graph G2(0);
    GraphUpdater gu2(G2);
    DynamicForestFireGenerator ffg2(1.0, true, 1.0);
    stream = ffg2.generate(10);
    gu2.update(stream);
    EXPECT_TRUE(G2.checkConsistency());
    EXPECT_EQ(10u, G2.numberOfNodes());
    G2.forNodePairs([&](node u, node v) {
        if (v < u) {
            EXPECT_TRUE(G2.hasEdge(u,v));
        }
    });
    stream = ffg2.generate(10);
    gu2.update(stream);
    EXPECT_EQ(20u, G2.numberOfNodes());
}

TEST_F(GeneratorsGTest, testRegularRingLatticeGenerator) {
    int n0 = 10;
    int neighbors = 2;
    auto testRingLattice = [&](Graph G) {
        EXPECT_EQ(n0, (int) G.numberOfNodes());
        EXPECT_EQ(n0 * neighbors, (int) G.numberOfEdges());
        G.forNodePairs([&](node u, node v) {
            int diff = std::abs((int) u- (int) v);
            if (u != v && (diff <= neighbors || diff >= n0 - neighbors)) {
                EXPECT_TRUE(G.hasEdge(u,v));
            } else {
                EXPECT_FALSE(G.hasEdge(u,v));
            }
        });
    };

    RegularRingLatticeGenerator rrlg = RegularRingLatticeGenerator(n0, neighbors);
    testRingLattice(rrlg.generate());
}

TEST_F(GeneratorsGTest, testWattsStrogatzGenerator) {
    int n0 = 10;
    int neighbors = 2;
    auto testRingLattice = [&](Graph G) {
        G.forNodePairs([&](node u, node v) {
            int diff = std::abs((int) u- (int) v);
            if (u != v && (diff <= neighbors || diff >= n0 - neighbors)) {
                EXPECT_TRUE(G.hasEdge(u,v));
            } else {
                EXPECT_FALSE(G.hasEdge(u,v));
            }
        });
    };

    WattsStrogatzGenerator wsg1 = WattsStrogatzGenerator(n0, neighbors, 0.0);
    testRingLattice(wsg1.generate());

    WattsStrogatzGenerator wsg2 = WattsStrogatzGenerator(n0, neighbors, 0.3);
    Graph G = wsg2.generate();
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_EQ(n0, (int) G.numberOfNodes());
    EXPECT_EQ(n0*neighbors, (int) G.numberOfEdges());
}

TEST_F(GeneratorsGTest, testWattsStrogatzGeneratorBigKs) {
    constexpr count nodes = 10;
    constexpr count neighbors = 4;
    constexpr double p = 0.5;
    for (int seed : {1, 2, 3, 4, 5}) {
        Aux::Random::setSeed(seed, false);
        const auto G = WattsStrogatzGenerator(nodes, neighbors, p).generate();
        EXPECT_TRUE(G.checkConsistency());
        EXPECT_EQ(nodes, G.numberOfNodes());
        EXPECT_EQ(nodes*neighbors, G.numberOfEdges());
    }

}

TEST_F(GeneratorsGTest, testDorogovtsevMendesGenerator) {
    int n0 = 20;
    DorogovtsevMendesGenerator dmg = DorogovtsevMendesGenerator(n0);
    Graph G = dmg.generate();

    EXPECT_EQ(n0, (int) G.numberOfNodes());
    EXPECT_EQ(2 * n0 - 3, (int) G.numberOfEdges());
    G.forNodes([&](node u) {
        count c = 0;
        G.forNeighborsOf(u, [&](node v) {
            if (v < u) {
                c += 1;
            }
        });
        if (u <= 2) {
            EXPECT_EQ(u, c);
        } else {
            EXPECT_EQ(2u, c);
        }
    });
    EXPECT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testDynamicDorogovtsevMendesGenerator) {
    count n0 = 20;
    DynamicDorogovtsevMendesGenerator ddmg = DynamicDorogovtsevMendesGenerator();
    Graph G(0);
    GraphUpdater gu(G);
    std::vector<GraphEvent> stream;
    stream = ddmg.generate(n0 - 3);
    gu.update(stream);

    EXPECT_EQ(n0, G.numberOfNodes());
    EXPECT_EQ(2*n0-3, G.numberOfEdges());
    G.forNodes([&](node u) {
        count c = 0;
        G.forNeighborsOf(u, [&](node v) {
            if (v < u) {
                c += 1;
            }
        });
        if (u <= 2) {
            EXPECT_EQ(u, c);
        } else {
            EXPECT_EQ(2u, c);
        }
    });
}

TEST_F(GeneratorsGTest, testStaticDegreeSequenceGenerator) {
    auto test_known = [](const std::vector<count> &seq, bool result) {
        HavelHakimiGenerator gen(seq, true);
        ASSERT_EQ(gen.isRealizable(), result);
    };

    // some contrived sequence
    test_known({1}, false);
    test_known({1, 1}, true);
    test_known({2, 2, 2}, true);
    test_known({1, 2, 1}, true);
    test_known({1, 3, 1}, false);
    test_known({1, 1, 3, 1}, true);

    // some random sequences
    std::mt19937_64 prng(1);
    std::uniform_int_distribution<node> distr_num_nodes(1, 50);

    unsigned numRealized = 0;
    for (int iter = 0; iter < 100; ++iter) {
        const auto n = distr_num_nodes(prng);
        std::vector<count> seq(n);
        std::generate(seq.begin(), seq.end(), [&] {
            return std::uniform_int_distribution<count>{0, n - 1}(prng);
        });

        HavelHakimiGenerator gen(seq, true);
        const auto isRealizable = gen.isRealizable();

        const auto G = gen.generate();
        const auto nodes = G.nodeRange();
        ASSERT_EQ(G.numberOfNodes(), n);

        const auto didRealize =
            std::all_of(nodes.begin(), nodes.end(), [&](node u) { return G.degree(u) == seq[u]; });

        ASSERT_EQ(isRealizable, didRealize);
        numRealized += didRealize;
    }

    ASSERT_GT(numRealized, 10);
    ASSERT_LT(numRealized, 90);
}

TEST_F(GeneratorsGTest, testStochasticBlockmodel) {
    count n = 10;
    count nBlocks = 2;
    std::vector<index> membership = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
    std::vector<std::vector<double> > affinity = {{1.0, 0.0}, {0.0, 1.0}};
    StochasticBlockmodel sbm(n, nBlocks, membership, affinity);
    Graph G = sbm.generate();

    EXPECT_EQ(n, G.numberOfNodes());
    EXPECT_EQ(20u, G.numberOfEdges());
}

/**
 * Test whether points generated in hyperbolic space fulfill basic constraints
 */
TEST_F(GeneratorsGTest, testHyperbolicPointGeneration) {
    count n = 1000;
    double stretch = Aux::Random::real(0.5,1.5);
    double alpha = Aux::Random::real(0.5,1.5);
    double R = HyperbolicSpace::hyperbolicAreaToRadius(n)*stretch;
    vector<double> angles(n, -1);
    vector<double> radii(n, -1);
    HyperbolicSpace::fillPoints(angles, radii, R, alpha);
    for (index i = 0; i < n; i++) {
        EXPECT_GE(angles[i], 0);
        EXPECT_LT(angles[i], 2*PI);
        EXPECT_GE(radii[i], 0);
        EXPECT_LE(radii[i], R);
    }
}

/**
 * Test whether the number edges generated by the hyperbolic generator agree at least roughly with theory
 */
TEST_F(GeneratorsGTest, testHyperbolicGenerator) {
    Aux::Random::setSeed(0, false);
    count n = 5000;
    double k = 16;
    count m = k*n/2;
    HyperbolicGenerator gen(n,k,7);
    Graph G = gen.generate();
    EXPECT_EQ(G.numberOfNodes(), n);
    EXPECT_TRUE(G.checkConsistency());
    EXPECT_NEAR(G.numberOfEdges(), m, m/5);
}

/**
 * Check consistency of graphs generated by the hyperbolic generator
 */
TEST_F(GeneratorsGTest, testHyperbolicGeneratorConsistency) {
    Aux::Random::setSeed(0, false);
    count n = 5000;
    double k = 6;
    count m = n*k/2;
    HyperbolicGenerator gen(n, k);
    Graph G = gen.generate();
    EXPECT_NEAR(G.numberOfEdges(), m, m/5);
    ASSERT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testHyperbolicGeneratorMechanicGraphs) {
    Aux::Random::setSeed(0, false);
    count n = 2000;
    double k = 6;
    count m = n*k/2;
    HyperbolicGenerator gen(n, k, 3, 0.14);
    Graph G = gen.generate();
    EXPECT_NEAR(G.numberOfEdges(), m, m/10);
    ASSERT_TRUE(G.checkConsistency());
}

TEST_F(GeneratorsGTest, testConfigurationModelGeneratorOnRealSequence) {
    METISGraphReader reader;
    std::vector<std::string> graphs = {"input/jazz.graph",
            "input/lesmis.graph"}; //, "input/PGPgiantcompo.graph", "input/coAuthorsDBLP.graph"};

    for (auto path : graphs) {
        Graph G = reader.read(path);
        count n = G.numberOfNodes();
        std::vector<count> sequence(n);
        G.forNodes([&](node u){
            sequence[u] = G.degree(u);
        });

        bool skipTest = false;
        EdgeSwitchingMarkovChainGenerator gen(sequence, skipTest);
        Graph G2 = gen.generate();

        count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
        EXPECT_EQ(volume, 2 * G2.numberOfEdges());

        if (volume < 50000) {
            std::vector<count> testSequence(n);
            G2.forNodes([&](node u){
                testSequence[u] = G2.degree(u);
            });
            Aux::Parallel::sort(testSequence.begin(), testSequence.end(), std::greater<count>());
            Aux::Parallel::sort(sequence.begin(), sequence.end(), std::greater<count>());

            for (index i = 0; i < n; ++i) {
                EXPECT_EQ(sequence[i], testSequence[i]);
            }
        }
    }
}

TEST_F(GeneratorsGTest, debugHyperbolicHighTemperatureGraphs) {
    count n = 10000;
    double k = 10;
    double gamma = 3;
    count m = n*k/2;
    for (double T = 0; T < 10; T += 0.1) {
        if (std::abs(T-1) < 0.00001) continue;
        HyperbolicGenerator gen(n, k, gamma, T);
        Graph G = gen.generate();
        EXPECT_NEAR(G.numberOfEdges(), m, m/10);
    }
}

TEST_F(GeneratorsGTest, debugGiganticCollectionOfHyperbolicTemperatureGraphs) {
    for (index i = 0; i < 30; i++) {
        count n = 10000;
        double k = 10;
        double T = 0.1;
        count m = n*k/2;
        HyperbolicGenerator gen(n, k, 3, T);
        Graph G = gen.generate();
        EXPECT_NEAR(G.numberOfEdges(), m, m/10);
        //EXPECT_TRUE(G.checkConsistency());
    }
}

TEST_F(GeneratorsGTest, debugGiganticCollectionOfHyperbolicUnitDiskGraphs) {
    count n = 1000000;
    double k = 1;
    for (index i = 0; i < 7; i++) {
        count m = n*k/2;
        HyperbolicGenerator gen(n, k, 7);
        Graph G = gen.generate();
        EXPECT_NEAR(G.numberOfEdges(), m, m/5);
        EXPECT_TRUE(G.checkConsistency());
        k *= 2;
    }
}

TEST_F(GeneratorsGTest, testLFRGenerator) {
    Aux::Random::setSeed(42, false);
    count n = 500;
    LFRGenerator gen(n);
    gen.generatePowerlawDegreeSequence(20, 50, -2);
    gen.generatePowerlawCommunitySizeSequence(10, 50, -1);
    gen.setMu(0.5);
    gen.run();
    Graph G1 = gen.getMoveGraph();
    gen.run(); // should rewire the edges but nothing else
    Graph G2 = gen.getMoveGraph();
    EXPECT_EQ(n, G1.numberOfNodes());
    EXPECT_EQ(n, G2.numberOfNodes());
    EXPECT_EQ(G1.numberOfEdges(), G2.numberOfEdges());
}

TEST_F(GeneratorsGTest, testLFRGeneratorImpossibleSequence) {
    LFRGenerator gen(100);
    gen.generatePowerlawDegreeSequence(10, 11, -2);
    EXPECT_ANY_THROW(gen.generatePowerlawCommunitySizeSequence(9, 8, -3));
    gen.setMu(0.5);
    EXPECT_THROW(gen.run(), std::runtime_error);
    EXPECT_THROW(gen.getMoveGraph(), std::runtime_error);
}

TEST_F(GeneratorsGTest, testLFRGeneratorWithRealData) {
    std::vector<count> degreeSequence = {124, 54, 7, 15, 81, 35, 4, 37, 3, 33, 2, 79, 2, 47, 40, 102, 17, 75, 2, 43, 49, 11, 12, 3, 87, 44, 67, 18, 22, 51, 96, 9, 17, 1, 1, 81, 3,
    44, 59, 2, 30, 85, 69, 28, 45, 12, 30, 38, 32, 20, 11, 28, 42, 30, 8, 86, 57, 56, 50, 51, 39, 80, 47, 124, 41, 15, 5, 3, 100, 57, 34, 37, 7, 2, 1, 11, 75, 137, 47, 11, 67, 92,
    67, 34, 156, 33, 8, 29, 35, 203, 20, 2, 51, 2, 54, 69, 14, 4, 1, 1, 72, 22, 45, 11, 34, 2, 21, 7, 4, 25, 1, 16, 63, 86, 25, 143, 45, 49, 23, 80, 3, 21, 1, 11, 14, 38, 73,
    14, 9, 45, 45, 117, 58, 5, 10, 33, 25, 28, 20, 41, 37, 95, 15, 57, 114, 84, 1, 22, 61, 22, 126, 132, 40, 77, 84, 64, 12, 95, 2, 39, 67, 40, 1, 30, 80, 57, 62, 48, 1, 64, 56,
    3, 26, 46, 79, 53, 38, 16, 26, 71, 3, 1, 49, 1, 18, 62, 39, 117, 9, 81, 50, 38, 4, 9, 68, 76, 61, 51, 49, 50, 8, 84, 56, 2, 59, 2, 1, 29, 28, 109, 33, 12, 37, 45, 12, 39,
    26, 42, 18, 54, 1, 11, 194, 39, 24, 65, 30, 59, 56, 16, 16, 41, 12, 30, 26, 60, 4, 13, 43, 14, 44, 62, 63, 60, 54, 9, 27, 53, 85, 33, 98, 107, 54, 7, 106, 38, 172, 18,
    30, 32, 56, 2, 58, 1, 132, 56, 30, 61, 2, 54, 31, 1, 13, 9, 113, 83, 27, 171, 120, 74, 25, 22, 48, 24, 51, 46, 26, 9, 31, 15, 58, 12, 33, 39, 56, 77, 1, 51, 11, 2, 89, 40,
    15, 47, 52, 1, 39, 104, 41, 45, 14, 57, 17, 53, 45, 10, 1, 1, 36, 40, 79, 45, 76, 3, 64, 3, 6, 45, 6, 153, 24, 47, 62, 52, 80, 23, 4, 51, 108, 4, 70, 80, 40, 37, 38, 20, 7,
    30, 7, 5, 9, 59, 35, 27, 144, 22, 23, 10, 9, 1, 98, 19, 2, 53, 37, 41, 53, 47, 22, 6, 63, 1, 142, 8, 99, 48, 144, 62, 28, 20, 67, 7, 84, 52, 28, 16, 8, 10, 65, 88, 26, 24, 49,
    1, 49, 48, 10, 2, 20, 102, 1, 3, 3, 63, 21, 42, 34, 21, 22, 2, 12, 22, 74, 9, 104, 53, 24, 68, 28, 132, 88, 40, 46, 138, 41, 2, 24, 2, 79, 11, 75, 46, 61, 59, 21, 92, 1, 32, 16,
    68, 9, 48, 66, 100, 58, 51, 35, 52, 14, 22, 23, 39, 121, 55, 105, 2, 38, 41, 70, 41, 13, 30, 110, 13, 82, 12, 25, 29, 59, 9, 6, 35, 63, 2, 20, 56, 97, 22, 44, 27, 135, 66, 85, 115,
    31, 40, 119, 24, 28, 65, 24, 9, 18, 103, 82, 13, 36, 102, 67, 41, 2, 104, 79, 4, 2, 11, 1, 35, 16, 28, 90, 62, 6, 119, 64, 100, 47, 20, 80, 55, 32, 45, 5, 89, 61, 75, 2, 58, 28, 35,
    1, 150, 1, 82, 34, 1, 65, 97, 61, 2, 44, 3, 35, 76, 6, 13, 7, 9, 67, 81, 6, 56, 5, 70, 70, 48, 50, 37, 96, 48, 47, 71, 138, 18, 31, 2, 88, 18, 15, 108, 65, 9, 9, 79, 4, 70, 14, 86,
    5, 16, 1, 44, 1, 32, 26, 22, 44, 23, 103, 104, 82, 115, 30, 9, 38, 43, 3, 38, 11, 26, 62, 40, 32, 29, 66, 9, 39, 4, 66, 14, 26, 19, 14, 33, 121, 13, 23, 20, 41, 40, 68, 10, 160, 60,
    105, 156, 4, 27, 2, 120, 136, 61, 2, 32, 11, 46, 23, 157, 20, 41, 64, 35, 23, 32, 36, 46, 79, 12, 68, 30, 68, 45, 34, 102, 20, 66, 76, 32, 76, 31, 7, 2, 13, 184, 33, 73, 18, 43, 35,
    141, 43, 77, 7, 9, 28, 39, 27, 44, 109, 4, 8, 41, 30, 3, 48, 8, 160, 33, 65, 14, 100, 34, 2, 94, 35, 19, 45, 19, 30, 59, 30, 11, 5, 1, 73, 19, 48, 1, 248, 15, 13, 3, 1, 36, 99, 95, 50,
    72, 57, 68, 11, 42, 8, 40, 5, 68, 12, 34, 15, 93, 14, 63, 22, 6, 152, 55, 35, 60, 14, 13, 5, 49, 2, 1, 79, 59, 29, 86, 5, 71, 2, 23, 9, 64, 33, 77, 19, 23, 61, 39, 9, 30, 55,
    23, 42, 16, 108, 42, 40};
    std::vector<count> partition = {1, 7, 3, 7, 1, 5, 13, 3, 7, 3, 8, 1, 16, 3, 2, 2, 10, 6, 11, 6, 7, 6, 3, 14, 3, 8, 2, 7, 2, 1, 5, 1, 4, 1, 18, 1, 14, 3, 1, 10, 4, 5, 1, 8, 7, 4,
    4, 4, 2, 10, 14, 2, 5, 6, 1, 5, 5, 4, 3, 1, 3, 6, 6, 5, 2, 8, 7, 6, 5, 5, 1, 1, 3, 16, 8, 2, 2, 6, 1, 12, 1, 3, 3, 5, 4, 6, 10, 3, 8, 1, 4, 4, 6, 4, 2, 1, 11, 4, 7, 1, 8, 3, 6,
    10, 4, 16, 1, 7, 6, 10, 2, 4, 4, 4, 6, 4, 4, 1, 6, 1, 12, 7, 1, 7, 8, 8, 5, 13, 6, 3, 4, 1, 5, 8, 7, 3, 9, 4, 9, 1, 3, 4, 2, 5, 4, 1, 18, 8, 3, 3, 4, 5, 1, 1, 2, 5, 3, 4, 8, 1,
    7, 7, 7, 7, 1, 3, 6, 3, 19, 1, 8, 7, 5, 7, 5, 3, 2, 5, 6, 3, 7, 3, 5, 5, 7, 5, 4, 2, 11, 4, 2, 5, 3, 6, 6, 2, 9, 1, 3, 5, 7, 8, 2, 10, 9, 17, 6, 8, 3, 6, 5, 3, 7, 9, 2, 2, 1, 1,
    4, 5, 13, 1, 3, 1, 10, 7, 5, 3, 1, 2, 3, 1, 3, 9, 6, 3, 11, 8, 2, 1, 5, 6, 9, 2, 6, 6, 2, 1, 1, 5, 5, 7, 7, 7, 6, 6, 5, 3, 6, 2, 2, 3, 2, 12, 4, 1, 6, 1, 13, 6, 1, 11, 1, 13, 3,
    1, 8, 6, 1, 1, 3, 2, 4, 9, 5, 6, 2, 7, 5, 2, 6, 8, 4, 2, 2, 8, 6, 9, 6, 8, 3, 6, 7, 4, 6, 1, 2, 8, 7, 3, 1, 1, 6, 5, 4, 10, 6, 6, 6, 6, 4, 8, 5, 1, 4, 1, 12, 3, 11, 1, 3, 10, 4,
    5, 1, 3, 6, 5, 2, 3, 1, 5, 1, 1, 6, 1, 4, 6, 6, 8, 7, 3, 9, 1, 7, 1, 1, 6, 8, 3, 9, 3, 6, 3, 7, 4, 10, 4, 7, 15, 3, 3, 4, 5, 1, 1, 1, 4, 5, 2, 2, 6, 8, 5, 4, 4, 12, 10, 6, 2, 9,
    4, 2, 7, 8, 7, 4, 6, 10, 2, 15, 12, 6, 3, 3, 7, 1, 3, 10, 4, 13, 2, 1, 3, 6, 3, 8, 2, 2, 4, 2, 6, 3, 2, 6, 4, 10, 4, 1, 3, 3, 7, 1, 6, 7, 1, 19, 4, 2, 6, 4, 7, 5, 2, 7, 2, 2, 3,
    2, 8, 2, 6, 4, 1, 1, 10, 7, 4, 3, 1, 9, 6, 5, 4, 2, 2, 1, 7, 1, 7, 8, 2, 2, 12, 1, 5, 1, 1, 6, 8, 1, 2, 1, 4, 6, 4, 2, 6, 8, 2, 2, 7, 1, 6, 5, 3, 4, 1, 1, 5, 17, 3, 1, 9, 4, 8,
    7, 8, 7, 4, 1, 2, 11, 2, 6, 5, 4, 7, 4, 3, 6, 7, 4, 5, 1, 5, 12, 5, 1, 2, 1, 2, 4, 8, 1, 4, 3, 6, 5, 12, 3, 9, 8, 2, 11, 6, 4, 7, 5, 4, 11, 1, 4, 1, 6, 3, 6, 9, 1, 4, 8, 2, 5, 4,
    5, 7, 2, 1, 3, 4, 1, 4, 7, 5, 8, 5, 8, 5, 8, 8, 2, 1, 10, 6, 7, 1, 8, 2, 3, 2, 6, 6, 4, 11, 3, 7, 3, 8, 10, 2, 3, 1, 3, 6, 4, 1, 4, 6, 6, 8, 4, 5, 3, 2, 5, 3, 4, 4, 5, 6, 10, 9,
    3, 4, 1, 2, 2, 5, 12, 3, 5, 4, 3, 8, 5, 4, 6, 4, 10, 5, 1, 4, 5, 1, 6, 3, 5, 4, 4, 6, 5, 1, 8, 4, 3, 3, 3, 6, 4, 4, 1, 7, 4, 8, 4, 2, 3, 4, 1, 2, 5, 4, 11, 2, 7, 1, 2, 5, 3, 12,
    4, 3, 3, 7, 1, 7, 3, 4, 7, 7, 5, 1, 2, 7, 3, 2, 5, 4, 4, 3, 9, 1, 4, 1, 3, 6, 8, 6, 2, 9, 4, 7, 2, 6, 3, 3, 6, 3, 2, 4, 4, 2, 6, 2, 4, 5, 2, 4, 5, 6, 3, 1, 7, 4, 4, 2, 3, 1, 8,
    3, 7, 3, 15, 13, 1, 8, 2, 3, 2, 6, 14, 1, 4, 1, 5, 5, 2, 5, 3, 9, 8, 7, 5, 2, 2, 2, 3, 6, 1};
    std::vector<double> mu = {0.6209677419354839, 0.6851851851851851, 0.2857142857142857, 0.6, 0.5308641975308642, 0.2857142857142857, 0.75, 0.3783783783783784, 0.6666666666666667,
    0.4545454545454546, 0.5, 0.4177215189873418, 0.0, 0.3191489361702128, 0.275, 0.5196078431372548, 0.3529411764705882, 0.7066666666666667, 0.5, 0.3023255813953488, 0.7142857142857143,
    0.2727272727272727, 0.08333333333333337, 0.0, 0.3448275862068966, 0.5909090909090908, 0.5970149253731343, 0.5, 0.2727272727272727, 0.2941176470588235, 0.375, 0.0, 0.23529411764705888,
    0.0, 0.0, 0.308641975308642, 0.33333333333333337, 0.09090909090909094, 0.4067796610169492, 0.5, 0.16666666666666663, 0.3529411764705882, 0.21739130434782605, 0.3928571428571429,
    0.6444444444444444, 0.6666666666666667, 0.4, 0.07894736842105265, 0.28125, 0.30000000000000004, 0.7272727272727273, 0.4642857142857143, 0.26190476190476186, 0.5, 0.375,
    0.36046511627906974, 0.3508771929824561, 0.2678571428571429, 0.12, 0.196078431372549, 0.23076923076923073, 0.4125, 0.34042553191489366, 0.6532258064516129, 0.31707317073170727,
    0.5333333333333333, 0.4, 0.0, 0.62, 0.08771929824561409, 0.4117647058823529, 0.16216216216216217, 0.1428571428571429, 0.0, 0.0, 0.2727272727272727, 0.52, 0.5547445255474452,
    0.3191489361702128, 0.36363636363636365, 0.35820895522388063, 0.7282608695652174, 0.34328358208955223, 0.08823529411764708, 0.7371794871794872, 0.1515151515151515, 0.125,
    0.5862068965517242, 0.6, 0.5172413793103448, 0.19999999999999996, 0.5, 0.4901960784313726, 0.0, 0.7777777777777778, 0.2028985507246377, 0.7142857142857143, 0.75, 0.0, 0.0,
    0.75, 0.2727272727272727, 0.37777777777777777, 0.09090909090909094, 0.20588235294117652, 0.0, 0.4285714285714286, 0.5714285714285714, 0.75, 0.43999999999999995, 0.0, 0.625,
    0.23809523809523814, 0.5813953488372092, 0.19999999999999996, 0.4475524475524476, 0.19999999999999996, 0.12244897959183676, 0.5217391304347826, 0.13749999999999996, 0.33333333333333337,
    0.47619047619047616, 0.0, 0.36363636363636365, 0.3571428571428571, 0.4473684210526315, 0.1917808219178082, 0.7142857142857143, 0.2222222222222222, 0.5111111111111111, 0.3111111111111111,
    0.23931623931623935, 0.4655172413793104, 0.8, 0.8, 0.5757575757575757, 0.64, 0.1428571428571429, 0.8, 0.31707317073170727, 0.5135135135135135, 0.5789473684210527, 0.06666666666666665,
    0.5263157894736843, 0.7192982456140351, 0.40476190476190477, 0.0, 0.7727272727272727, 0.3114754098360656, 0.13636363636363635, 0.5714285714285714, 0.5681818181818181, 0.475,
    0.1428571428571429, 0.5357142857142857, 0.40625, 0.16666666666666663, 0.49473684210526314, 0.5, 0.3846153846153846, 0.5970149253731343, 0.625, 0.0, 0.5, 0.30000000000000004,
    0.4035087719298246, 0.3548387096774194, 0.10416666666666663, 0.0, 0.5625, 0.5178571428571428, 0.6666666666666667, 0.23076923076923073, 0.6521739130434783, 0.6075949367088608,
    0.7924528301886793, 0.3157894736842105, 0.375, 0.6923076923076923, 0.676056338028169, 0.0, 0.0, 0.44897959183673475, 0.0, 0.4444444444444444, 0.27419354838709675, 0.5641025641025641,
    0.3931623931623932, 0.5555555555555556, 0.49382716049382713, 0.14, 0.3157894736842105, 0.0, 0.4444444444444444, 0.47058823529411764, 0.5, 0.8524590163934427, 0.4117647058823529,
    0.22448979591836737, 0.16000000000000003, 0.5, 0.7261904761904762, 0.2857142857142857, 0.5, 0.8305084745762712, 0.5, 0.0, 0.48275862068965514, 0.3571428571428571, 0.5688073394495412,
    0.24242424242424243, 0.41666666666666663, 0.6486486486486487, 0.7333333333333334, 0.33333333333333337, 0.33333333333333337, 0.23076923076923073, 0.19047619047619047,
    0.16666666666666663, 0.40740740740740744, 0.0, 0.4545454545454546, 0.7371134020618557, 0.2564102564102564, 0.5416666666666667, 0.8, 0.3666666666666667, 0.30508474576271183,
    0.4464285714285714, 0.125, 0.3125, 0.24390243902439024, 0.25, 0.7, 0.11538461538461542, 0.4, 0.5, 0.46153846153846156, 0.39534883720930236, 0.2857142857142857, 0.31818181818181823,
    0.467741935483871, 0.8253968253968254, 0.01666666666666672, 0.6666666666666667, 0.4444444444444444, 0.4814814814814815, 0.41509433962264153, 0.6, 0.36363636363636365, 0.4387755102040817,
    0.6448598130841121, 0.537037037037037, 0.4285714285714286, 0.5377358490566038, 0.2894736842105263, 0.6046511627906976, 0.11111111111111116, 0.6666666666666667, 0.1875, 0.3928571428571429,
    0.5, 0.3448275862068966, 0.0, 0.5833333333333333, 0.5178571428571428, 0.30000000000000004, 0.180327868852459, 0.5, 0.42592592592592593, 0.32258064516129037, 0.0, 0.07692307692307687,
    0.8888888888888888, 0.6283185840707964, 0.5060240963855422, 0.5185185185185186, 0.6842105263157895, 0.3916666666666667, 0.6486486486486487, 0.43999999999999995, 0.13636363636363635,
    0.375, 0.7916666666666666, 0.2941176470588235, 0.5869565217391304, 0.23076923076923073, 0.4444444444444444, 0.3870967741935484, 0.4, 0.5, 0.41666666666666663, 0.5454545454545454,
    0.15384615384615385, 0.625, 0.7532467532467533, 0.0, 0.803921568627451, 0.36363636363636365, 0.0, 0.4606741573033708, 0.625, 0.19999999999999996, 0.3829787234042553, 0.6346153846153846,
    0.0, 0.2564102564102564, 0.6442307692307692, 0.8780487804878049, 0.37777777777777777, 0.0, 0.4035087719298246, 0.4117647058823529, 0.3584905660377359, 0.28888888888888886,
    0.19999999999999996, 0.0, 0.0, 0.25, 0.275, 0.49367088607594933, 0.7555555555555555, 0.5394736842105263, 0.0, 0.5625, 0.33333333333333337, 0.6666666666666667, 0.3555555555555555,
    0.5, 0.5163398692810457, 0.125, 0.7659574468085106, 0.33870967741935487, 0.23076923076923073, 0.3375, 0.5217391304347826, 0.5, 0.3137254901960784, 0.5555555555555556, 0.5,
    0.5571428571428572, 0.75, 0.25, 0.21621621621621623, 0.21052631578947367, 0.19999999999999996, 0.1428571428571429, 0.5333333333333333, 0.5714285714285714, 0.6, 0.2222222222222222,
    0.423728813559322, 0.7714285714285715, 0.14814814814814814, 0.7847222222222222, 0.2272727272727273, 0.21739130434782605, 0.19999999999999996, 0.5555555555555556, 0.0, 0.8775510204081632,
    0.3157894736842105, 0.0, 0.6226415094339622, 0.5405405405405406, 0.6341463414634146, 0.6981132075471699, 0.12765957446808507, 0.6818181818181819, 0.6666666666666667, 0.3492063492063492,
    0.0, 0.6619718309859155, 0.625, 0.4242424242424242, 0.22916666666666663, 0.5138888888888888, 0.5, 0.3214285714285714, 0.050000000000000044, 0.4626865671641791, 0.2857142857142857,
    0.7380952380952381, 0.6153846153846154, 0.3928571428571429, 0.1875, 0.5, 0.09999999999999998, 0.4769230769230769, 0.3522727272727273, 0.6923076923076923, 0.20833333333333337,
    0.24489795918367352, 0.0, 0.6326530612244898, 0.625, 0.0, 0.5, 0.5, 0.5686274509803921, 0.0, 0.33333333333333337, 0.6666666666666667, 0.31746031746031744, 0.19047619047619047,
    0.7380952380952381, 0.2647058823529411, 0.5714285714285714, 0.5, 0.0, 0.75, 0.09090909090909094, 0.3918918918918919, 0.33333333333333337, 0.5096153846153846, 0.13207547169811318,
    0.75, 0.4852941176470589, 0.1428571428571429, 0.49242424242424243, 0.43181818181818177, 0.275, 0.28260869565217395, 0.644927536231884, 0.2682926829268293, 0.0, 0.5, 0.5,
    0.430379746835443, 0.4545454545454546, 0.72, 0.4782608695652174, 0.4426229508196722, 0.35593220338983056, 0.38095238095238093, 0.44565217391304346, 0.0, 0.3125, 0.125,
    0.5735294117647058, 0.4444444444444444, 0.6458333333333333, 0.4545454545454546, 0.5700000000000001, 0.7758620689655172, 0.37254901960784315, 0.34285714285714286, 0.6153846153846154,
    0.5, 0.4545454545454546, 0.4782608695652174, 0.3076923076923077, 0.6611570247933884, 0.4545454545454546, 0.2952380952380952, 0.5, 0.736842105263158, 0.3902439024390244,
    0.6142857142857143, 0.41463414634146345, 0.6923076923076923, 0.5333333333333333, 0.40909090909090906, 0.07692307692307687, 0.3292682926829268, 0.5833333333333333, 0.6,
    0.4137931034482759, 0.3728813559322034, 0.4444444444444444, 0.33333333333333337, 0.3142857142857143, 0.6984126984126984, 0.0, 0.55, 0.5714285714285714, 0.4639175257731959,
    0.13636363636363635, 0.7045454545454546, 0.7037037037037037, 0.6370370370370371, 0.5151515151515151, 0.3529411764705882, 0.4782608695652174, 0.7096774193548387, 0.4, 0.7058823529411764,
    0.5416666666666667, 0.5357142857142857, 0.6153846153846154, 0.29166666666666663, 0.6666666666666667, 0.5555555555555556, 0.7184466019417476, 0.47560975609756095, 0.23076923076923073,
    0.13888888888888884, 0.4411764705882353, 0.25373134328358204, 0.1707317073170732, 0.5, 0.40384615384615385, 0.6075949367088608, 0.5, 0.5, 0.09090909090909094, 0.0, 0.6, 0.6875,
    0.1785714285714286, 0.5444444444444445, 0.33870967741935487, 0.33333333333333337, 0.4369747899159664, 0.28125, 0.64, 0.14893617021276595, 0.44999999999999996, 0.35, 0.6, 0.75,
    0.5111111111111111, 0.6, 0.6179775280898876, 0.4590163934426229, 0.29333333333333333, 0.0, 0.31034482758620685, 0.3928571428571429, 0.4, 0.0, 0.6266666666666667, 0.0, 0.7560975609756098,
    0.38235294117647056, 0.0, 0.3076923076923077, 0.6701030927835052, 0.2295081967213115, 0.5, 0.2272727272727273, 0.33333333333333337, 0.6857142857142857, 0.42105263157894735,
    0.33333333333333337, 0.5384615384615384, 0.5714285714285714, 0.4444444444444444, 0.34328358208955223, 0.4814814814814815, 0.16666666666666663, 0.5535714285714286, 0.19999999999999996,
    0.6428571428571428, 0.6857142857142857, 0.125, 0.21999999999999997, 0.8648648648648649, 0.42708333333333337, 0.5833333333333333, 0.5531914893617021, 0.3380281690140845,
    0.572463768115942, 0.38888888888888884, 0.5483870967741935, 0.0, 0.5681818181818181, 0.16666666666666663, 0.1333333333333333, 0.5555555555555556, 0.2153846153846154, 0.0,
    0.33333333333333337, 0.30379746835443033, 0.75, 0.6714285714285715, 0.5, 0.4534883720930233, 0.19999999999999996, 0.5, 0.0, 0.25, 0.0, 0.25, 0.2692307692307693, 0.18181818181818177,
    0.5, 0.30434782608695654, 0.4563106796116505, 0.29807692307692313, 0.4512195121951219, 0.6782608695652175, 0.30000000000000004, 0.6666666666666667, 0.42105263157894735,
    0.7674418604651163, 0.0, 0.736842105263158, 0.2727272727272727, 0.3076923076923077, 0.19354838709677424, 0.35, 0.125, 0.5517241379310345, 0.19696969696969702, 0.11111111111111116,
    0.4871794871794872, 0.5, 0.6515151515151515, 0.7857142857142857, 0.2692307692307693, 0.26315789473684215, 0.0714285714285714, 0.5757575757575757, 0.6859504132231404,
    0.23076923076923073, 0.04347826086956519, 0.35, 0.5609756097560976, 0.22499999999999998, 0.8529411764705882, 0.8, 0.775, 0.5, 0.4095238095238095, 0.6346153846153846, 0.25,
    0.33333333333333337, 0.0, 0.6416666666666666, 0.5220588235294117, 0.3770491803278688, 0.5, 0.75, 0.2727272727272727, 0.7173913043478262, 0.34782608695652173, 0.6178343949044586,
    0.4, 0.19512195121951215, 0.234375, 0.37142857142857144, 0.6086956521739131, 0.6875, 0.5555555555555556, 0.4130434782608695, 0.3417721518987342, 0.0, 0.5441176470588236,
    0.5, 0.4117647058823529, 0.4666666666666667, 0.4411764705882353, 0.48039215686274506, 0.30000000000000004, 0.6363636363636364, 0.6447368421052632, 0.65625, 0.5263157894736843,
    0.4516129032258065, 0.1428571428571429, 0.5, 0.23076923076923073, 0.8206521739130435, 0.696969696969697, 0.4794520547945206, 0.6666666666666667, 0.2325581395348837, 0.4571428571428572,
    0.6382978723404256, 0.2558139534883721, 0.6753246753246753, 0.5714285714285714, 0.2222222222222222, 0.6785714285714286, 0.33333333333333337, 0.14814814814814814, 0.6136363636363636,
    0.4128440366972477, 0.75, 0.5, 0.24390243902439024, 0.16666666666666663, 0.6666666666666667, 0.375, 0.75, 0.65625, 0.5757575757575757, 0.6153846153846154, 0.5, 0.64, 0.08823529411764708,
    0.0, 0.7765957446808511, 0.6857142857142857, 0.3157894736842105, 0.2666666666666667, 0.3157894736842105, 0.3666666666666667, 0.576271186440678, 0.7666666666666666, 0.18181818181818177,
    0.6, 0.0, 0.17808219178082196, 0.5789473684210527, 0.6666666666666667, 0.0, 0.6774193548387097, 0.6666666666666667, 0.6923076923076923, 0.6666666666666667, 0.0, 0.19444444444444442,
    0.5858585858585859, 0.5684210526315789, 0.36, 0.1527777777777778, 0.24561403508771928, 0.38235294117647056, 0.5454545454545454, 0.23809523809523814, 0.25, 0.4, 0.19999999999999996,
    0.16176470588235292, 0.5, 0.17647058823529416, 0.4666666666666667, 0.3763440860215054, 0.1428571428571429, 0.47619047619047616, 0.7272727272727273, 0.33333333333333337,
    0.618421052631579, 0.34545454545454546, 0.37142857142857144, 0.3833333333333333, 0.3571428571428571, 0.0, 0.4, 0.34693877551020413, 0.5, 0.0, 0.49367088607594933, 0.6610169491525424,
    0.1724137931034483, 0.2093023255813954, 0.6, 0.6056338028169015, 0.0, 0.26086956521739135, 0.2222222222222222, 0.25, 0.3939393939393939, 0.3246753246753247, 0.21052631578947367,
    0.4782608695652174, 0.2622950819672131, 0.7435897435897436, 0.4444444444444444, 0.6, 0.2727272727272727, 0.4782608695652174, 0.26190476190476186, 0.1875, 0.6018518518518519,
    0.3571428571428571, 0.050000000000000044};
    Partition C(partition.size());
    C.setUpperBound(20);
    for (node u = 0; u < partition.size(); ++u) {
        C[u] = partition[u];
    }
    LFRGenerator gen(degreeSequence.size());
    gen.setDegreeSequence(degreeSequence);
    gen.setPartition(C);
    gen.setMu(mu);
    gen.run();
    Graph G = gen.getGraph();
    G.parallelForNodes([&](node u){
        EXPECT_EQ(G.degree(u),degreeSequence[u]);
    });
    EXPECT_EQ(C.numberOfSubsets(),gen.getPartition().numberOfSubsets());
}

TEST_F(GeneratorsGTest, testMocnikGenerator) {
    count dim = 3;
    count n = 10000;
    double k = 2.6;

    MocnikGenerator Mocnik(dim, n, k);
    Graph G(0);
    EXPECT_TRUE(G.isEmpty());
    G = Mocnik.generate();
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(G.numberOfNodes(), n);
    EXPECT_NEAR(G.numberOfEdges() * 1. / G.numberOfNodes(), std::pow(k, dim), 20000);
}

TEST_F(GeneratorsGTest, testMocnikGeneratorBasic) {
    count dim = 3;
    count n = 5000;
    double k = 2.6;

    MocnikGenerator Mocnik(dim, n, k);
    Graph G(0);
    EXPECT_TRUE(G.isEmpty());
    G = Mocnik.generate();
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(G.numberOfNodes(), n);
    EXPECT_NEAR(G.numberOfEdges() * 1. / G.numberOfNodes(), std::pow(k, dim), 10000);
}

} /* namespace NetworKit */
