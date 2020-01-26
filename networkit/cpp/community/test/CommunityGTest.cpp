/*
 * CommunityGTest.cpp
 *
 *  Created on: 27.02.2014
 *      Author: cls
 */

#include <gtest/gtest.h>

#include <tlx/unused.hpp>
#include <networkit/community/PLP.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/community/ParallelAgglomerativeClusterer.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/community/EdgeCut.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/community/GraphClusteringTools.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/community/Coverage.hpp>
#include <networkit/community/JaccardMeasure.hpp>
#include <networkit/community/NodeStructuralRandMeasure.hpp>
#include <networkit/community/GraphStructuralRandMeasure.hpp>
#include <networkit/community/NMIDistance.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/community/SampledGraphStructuralRandMeasure.hpp>
#include <networkit/community/SampledNodeStructuralRandMeasure.hpp>
#include <networkit/community/PartitionIntersection.hpp>
#include <networkit/community/HubDominance.hpp>
#include <networkit/community/IntrapartitionDensity.hpp>
#include <networkit/community/PartitionFragmentation.hpp>
#include <networkit/generators/ClusteredRandomGraphGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/community/CoverF1Similarity.hpp>
#include <networkit/structures/Cover.hpp>
#include <networkit/community/LPPotts.hpp>

namespace NetworKit {

class CommunityGTest: public testing::Test{};

TEST_F(CommunityGTest, testLabelPropagationOnUniformGraph) {
    ErdosRenyiGenerator graphGen(100, 0.2);
    Graph G = graphGen.generate();

    PLP lp(G);
    lp.run();
    Partition zeta = lp.getPartition();

    EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta)) << "the resulting partition should be a proper clustering";

    Modularity modularity;
    double mod = modularity.getQuality(zeta, G);
    DEBUG("modularity produced by LabelPropagation: " , mod);
    EXPECT_GE(1.0, mod) << "valid modularity values are in [-0.5, 1]";
    EXPECT_LE(-0.5, mod) << "valid modularity values are in [-0.5, 1]";
}


TEST_F(CommunityGTest, testLabelPropagationOnClusteredGraph_ForNumberOfClusters) {
    int64_t n = 100;
    count k = 3; // number of clusters

    ClusteredRandomGraphGenerator graphGen(n, k, 1.0, 0.00);
    Graph G = graphGen.generate();

    PLP lp(G);
    lp.run();
    Partition zeta = lp.getPartition();

    Modularity modularity;
    DEBUG("modularity produced by LabelPropagation: " , modularity.getQuality(zeta, G));

    EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta)) << "the resulting partition should be a proper clustering";
    EXPECT_EQ(k, zeta.numberOfSubsets()) << " " << k << " clusters are easy to detect";
}



TEST_F(CommunityGTest, testLabelPropagationOnDisconnectedGraph) {
    count n = 100;
    count k = 2; // number of clusters
    ClusteredRandomGraphGenerator graphGen(n, k, 1.0, 0.00);
    Graph G = graphGen.generate();


    PLP lp(G);
    lp.run();
    Partition zeta = lp.getPartition();

    Modularity modularity;
    DEBUG("modularity produced by LabelPropagation: " , modularity.getQuality(zeta, G));

    EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta)) << "the resulting partition should be a proper clustering";
    EXPECT_EQ(k, zeta.numberOfSubsets()) << " " << k << " clusters are easy to detect"; //FIXME

}


TEST_F(CommunityGTest, testLabelPropagationOnSingleNodeWithSelfLoop) {
    Graph G(1, true);
    node v = 0;
    G.setWeight(v, v, 42.0);

    PLP lp(G);
    lp.run();
    Partition zeta = lp.getPartition();

    EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta));
    EXPECT_TRUE(GraphClusteringTools::isSingletonClustering(G, zeta));
    EXPECT_TRUE(GraphClusteringTools::isOneClustering(G, zeta)); //FIXME does this make sense? singleton and one partition at the same time.

    Modularity modularity;
    DEBUG("modularity produced by LabelPropagation: " , modularity.getQuality(zeta, G));
}


TEST_F(CommunityGTest, testLabelPropagationOnManySmallClusters) {
    int64_t n = 1000;
    int k = 100; // number of clusters
    double pin = 1.0;
    double pout = 0.0;

    ClusteredRandomGraphGenerator graphGen(n, k, pin, pout);
  Aux::Random::setSeed(42, false);
    Graph G = graphGen.generate();
    Partition reference = graphGen.getCommunities();

    PLP lp(G);
    lp.run();
    Partition zeta = lp.getPartition();

    Modularity modularity;
    DEBUG("modularity produced by LabelPropagation: " , modularity.getQuality(zeta, G));
    DEBUG("number of clusters produced by LabelPropagation: k=" , zeta.numberOfSubsets());

    EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta)) << "the resulting partition should be a proper clustering";
    EXPECT_TRUE(GraphClusteringTools::equalClusterings(zeta, reference, G)) << "Can LabelPropagation detect the reference clustering?";

}



/*
TEST_F(CommunityGTest, testLouvainParallel2Naive) {
    count n = 1000;
    count k = 100;
    double pin = 1.0;
    double pout = 0.005;
    GraphGenerator graphGen;
    Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

    LouvainParallel louvain(G);
    louvain.run();
    Partition zeta = louvain.getPartition();

    INFO("number of clusters: " , zeta.numberOfSubsets());

    Modularity modularity;
    INFO("modularity: " , modularity.getQuality(zeta, G));

}
*/




TEST_F(CommunityGTest, testPLM) {
    METISGraphReader reader;
    Modularity modularity;
    Graph G = reader.read("input/PGPgiantcompo.graph");

    PLM plm(G, false, 1.0);
    plm.run();
    Partition zeta = plm.getPartition();

    DEBUG("number of clusters: " , zeta.numberOfSubsets());
    DEBUG("modularity: " , modularity.getQuality(zeta, G));
    EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta));

    PLM plmr(G, true, 1.0);
    plmr.run();
    Partition zeta2 = plmr.getPartition();

    DEBUG("number of clusters: " , zeta2.numberOfSubsets());
    DEBUG("modularity: " , modularity.getQuality(zeta2, G));
    EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta2));

}

TEST_F(CommunityGTest, testDeletedNodesPLM) {
    METISGraphReader reader;
    Modularity modularity;
    Graph G = reader.read("input/PGPgiantcompo.graph");

    G.forNeighborsOf(10, [&](node v) {
        G.removeEdge(10, v);
    });

    G.removeNode(10);

    PLM plm(G, false, 1.0);
    plm.run();
    Partition zeta = plm.getPartition();

    DEBUG("number of clusters: " , zeta.numberOfSubsets());
    DEBUG("modularity: " , modularity.getQuality(zeta, G));
    EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta));

    PLM plmr(G, true, 1.0);
    plmr.run();
    Partition zeta2 = plmr.getPartition();

    DEBUG("number of clusters: " , zeta2.numberOfSubsets());
    DEBUG("modularity: " , modularity.getQuality(zeta2, G));
    EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta2));

}

TEST_F(CommunityGTest, testModularity) {

    count n = 100;

    DEBUG("testing modularity on clustering of complete graph with " , n , " nodes");

    // make complete graph
    Graph G = Graph(n);
    G.forNodePairs([&](node u, node v) {
        G.addEdge(u,v);
    });
    DEBUG("total edge weight: " , G.totalEdgeWeight());

    ClusteringGenerator clusteringGenerator;

    Partition singleton = clusteringGenerator.makeSingletonClustering(G);
    Partition one = clusteringGenerator.makeOneClustering(G);

    Modularity modularity;

    DEBUG("calculating modularity for singleton clustering");
    double modSingleton = modularity.getQuality(singleton, G);

    DEBUG("calculating modularity for 1-clustering");
    double modOne = modularity.getQuality(one, G);

    DEBUG("mod(singleton-clustering) = " , modSingleton);
    DEBUG("mod(1-clustering) = " , modOne);


    EXPECT_EQ(0.0, modOne) << "1-clustering should have modularity of 0.0";
    EXPECT_GE(0.0, modSingleton) << "singleton clustering should have modularity less than 0.0";

}

TEST_F(CommunityGTest, testCoverage) {

    count n = 100;

    DEBUG("testing coverage on clustering of complete graph with " , n , " nodes");

    // make complete graph
    Graph G = Graph(n);
    G.forNodePairs([&](node u, node v) {
        G.addEdge(u,v);
    });

    ClusteringGenerator clusteringGenerator;

    Partition singleton = clusteringGenerator.makeSingletonClustering(G);
    Partition one = clusteringGenerator.makeOneClustering(G);

    Coverage coverage;

    DEBUG("calculating coverage for singleton clustering");
    double covSingleton = coverage.getQuality(singleton, G);

    DEBUG("calculating coverage for 1-clustering");
    double covOne = coverage.getQuality(one, G);

    DEBUG("mod(singleton-clustering) = " , covSingleton);
    DEBUG("mod(1-clustering) = " , covOne);


    EXPECT_EQ(1.0, covOne) << "1-clustering should have coverage of 1.0";
    EXPECT_GE(0.0, covSingleton) << "singleton clustering should have coverage of 0.0";

}

// TODO necessary testcase? move equals to some class ?
TEST_F(CommunityGTest, testClusteringEquality) {
    count n = 100;

    // make complete graph
    Graph G = Graph(n);
    G.forNodePairs([&](node u, node v) {
        G.addEdge(u,v);
    });

    ClusteringGenerator clusteringGen;
    Partition one1 = clusteringGen.makeOneClustering(G);
    Partition one2 = clusteringGen.makeOneClustering(G);

    EXPECT_TRUE(GraphClusteringTools::equalClusterings(one1, one2, G)) << "two 1-clusterings of G should be equal";

    Partition singleton1 = clusteringGen.makeSingletonClustering(G);
    Partition singleton2 = clusteringGen.makeSingletonClustering(G);

    EXPECT_TRUE(GraphClusteringTools::equalClusterings(singleton1, singleton2, G)) << "two singleton clusterings of G should be equal";

}


TEST_F(CommunityGTest, testEdgeCutMeasure) {
    /* Graph:
        0    3
         \  / \
          2    5
         /  \ /
        1    4
     */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);

    Partition part(n);
    part[0] = 0;
    part[1] = 0;
    part[2] = 0;
    part[3] = 1;
    part[4] = 2;
    part[5] = 1;

    EdgeCut ec;
    edgeweight cut = ec.getQuality(part, G);
    EXPECT_EQ(cut, 3);
}


TEST_F(CommunityGTest, testJaccardMeasure) {
    count n = 100;
    // make complete graph
    Graph G = Graph(n);
    G.forNodePairs([&](node u, node v) {
        G.addEdge(u,v);
    });

    ClusteringGenerator clusteringGen;
    Partition singleton = clusteringGen.makeSingletonClustering(G);
    Partition random = clusteringGen.makeRandomClustering(G, 10);

    JaccardMeasure jaccard;
    double j = jaccard.getDissimilarity(G, singleton, random);

    EXPECT_EQ(1.0, j) << "The singleton clustering and any other clustering compare with a dissimilarity of 1.0";

}


TEST_F(CommunityGTest, testNodeStructuralRandMeasure) {
    count n = 100;
    // make complete graph
    Graph G = Graph(n);
    G.forNodePairs([&](node u, node v) {
        G.addEdge(u,v);
    });

    ClusteringGenerator clusteringGen;
    Partition one1 = clusteringGen.makeOneClustering(G);
    Partition one2 = clusteringGen.makeOneClustering(G);

    NodeStructuralRandMeasure rand;
    double r = rand.getDissimilarity(G, one1, one2);

    EXPECT_EQ(0.0, r) << "Identical clusterings should compare with a dissimilarity of 0.0";

}

TEST_F(CommunityGTest, testGraphStructuralRandMeasure) {
    count n = 100;
    // make complete graph
    Graph G = Graph(n);
    G.forNodePairs([&](node u, node v) {
        G.addEdge(u,v);
    });

    ClusteringGenerator clusteringGen;
    Partition one1 = clusteringGen.makeOneClustering(G);
    Partition one2 = clusteringGen.makeOneClustering(G);

    GraphStructuralRandMeasure rand;
    double r = rand.getDissimilarity(G, one1, one2);

    EXPECT_EQ(0.0, r) << "Identical clusterings should compare with a dissimilarity of 0.0";

}

TEST_F(CommunityGTest, testNMIDistance) {
    // two 1-clusterings should have NMID = 0 because H is 0
    Graph G(1000);

    ClusteringGenerator clustGen;
    Partition one1 = clustGen.makeOneClustering(G);
    Partition one2 = clustGen.makeOneClustering(G);

    NMIDistance NMID;
    double distOne = NMID.getDissimilarity(G, one1, one2);

    DEBUG("NMID for two 1-clusterings: " , distOne);
    EXPECT_TRUE(Aux::NumericTools::equal(0.0, distOne)) << "NMID of two 1-clusterings should be 0.0";


    Partition singleton1 = clustGen.makeSingletonClustering(G);
    Partition singleton2 = clustGen.makeSingletonClustering(G);

    double distSingleton = NMID.getDissimilarity(G, singleton1, singleton2);
    DEBUG("NMID for two singleton clusterings: " , distSingleton);


    EXPECT_TRUE(Aux::NumericTools::equal(0.0, distSingleton)) << "NMID of two identical singleton clusterings should be 0.0";

    Partition continuous1 = clustGen.makeContinuousBalancedClustering(G, 40);
    Partition continuous2 = clustGen.makeContinuousBalancedClustering(G, 70);

    double distContinuous = NMID.getDissimilarity(G, continuous1, continuous2);
    DEBUG("NMID for two continuous clusterings: " , distContinuous);
    tlx::unused(distContinuous);

    Partition smallClusters = clustGen.makeContinuousBalancedClustering(G, 300);
    double distSingleIntersection = NMID.getDissimilarity(G, singleton1, smallClusters);
    EXPECT_LE(0.0, distSingleIntersection) << "NMID always needs to be 0 or positive";
}


TEST_F(CommunityGTest, testSampledRandMeasures) {
    count n = 42;
    // make complete graph
    Graph G = Graph(n);
    G.forNodePairs([&](node u, node v) {
        G.addEdge(u,v);
    });
    ClusteringGenerator clusteringGenerator;
    Partition one = clusteringGenerator.makeOneClustering(G);
    Partition singleton = clusteringGenerator.makeSingletonClustering(G);
    SampledNodeStructuralRandMeasure nRand(20);
    SampledGraphStructuralRandMeasure gRand(20);
    ASSERT_EQ(nRand.getDissimilarity(G, one, singleton),gRand.getDissimilarity(G, one, singleton));
}


TEST_F(CommunityGTest, debugParallelAgglomerativeAndPLM) {
    METISGraphReader reader;
    Graph jazz = reader.read("input/jazz.graph");
    Graph blog = reader.read("input/polblogs.graph");
    // FIXME: implementation of ParallelAgglomerativeClusterer needs overhaul
    Modularity modularity;
    ParallelAgglomerativeClusterer aggl(jazz);
    PLM louvain(jazz);
    // *** jazz graph
    // parallel agglomerative
    aggl.run();
    Partition clustering = aggl.getPartition();
    INFO("Match-AGGL number of jazz clusters: " , clustering.numberOfSubsets());
    INFO("Match-AGGL modularity jazz graph:   " , modularity.getQuality(clustering, jazz));

    // Louvain
    louvain.run();
    clustering = louvain.getPartition();
    INFO("Louvain number of jazz clusters: " , clustering.numberOfSubsets());
    INFO("Louvain modularity jazz graph:   " , modularity.getQuality(clustering, jazz));



    // *** blog graph
    ParallelAgglomerativeClusterer aggl2(jazz);
    PLM louvain2(jazz);
    // parallel agglomerative
    aggl2.run();
    clustering = aggl2.getPartition();
    INFO("Match-AGGL number of blog clusters: " , clustering.numberOfSubsets());
    INFO("Match-AGGL modularity blog graph:   " , modularity.getQuality(clustering, blog));

    // Louvain
    louvain2.run();
    clustering = louvain2.getPartition();
    INFO("Louvain number of blog clusters: " , clustering.numberOfSubsets());
    INFO("Louvain modularity blog graph:   " , modularity.getQuality(clustering, blog));
}

TEST_F(CommunityGTest, testClusteringIntersection) {
    PartitionIntersection intersection;
    // make complete graph
    count n = 1200;
    Graph G = Graph(n);
    G.forNodePairs([&](node u, node v) {
        G.addEdge(u,v);
    });

    ClusteringGenerator clusteringGenerator;
    Partition twelve = clusteringGenerator.makeContinuousBalancedClustering(G, 12);
    Partition singleton = clusteringGenerator.makeSingletonClustering(G);
    Partition eight = clusteringGenerator.makeContinuousBalancedClustering(G, 8);
    EXPECT_TRUE(GraphClusteringTools::equalClusterings(twelve, intersection.calculate(twelve, twelve), G)) << "Intersection of itself does not modify the clustering";
    EXPECT_TRUE(GraphClusteringTools::equalClusterings(singleton, intersection.calculate(twelve, singleton), G)) << "Intersection of singleton with any clustering is the singleton clustering";
    Partition sixteen = intersection.calculate(twelve, eight);
    EXPECT_EQ(16u, sixteen.numberOfSubsets());
    auto clusterSizes = sixteen.subsetSizeMap();
    size_t i = 0;
    for (auto size : clusterSizes) {
        if (i % 4 == 0 || i % 4 == 3) {
            EXPECT_EQ(100u, size.second) << "cluster size pattern is not 100 | 50 | 50 | 100 | 100 | 50 | 50 | 100 | 100 | ... | 50 | 50 | 100";
        } else {
            EXPECT_EQ(50u, size.second) << "cluster size pattern is not 100 | 50 | 50 | 100 | 100 | 50 | 50 | 100 | 100 | ... | 50 | 50 | 100";
        }
        ++i;
    }
}

TEST_F(CommunityGTest, testMakeNoncontinuousClustering) {
    ClusteringGenerator generator;
    // make complete graph
    count n = 100;
    Graph G = Graph(n);
    G.forNodePairs([&](node u, node v) {
        G.addEdge(u,v);
    });

    Partition con = generator.makeContinuousBalancedClustering(G, 10);
    Partition nonCon = generator.makeNoncontinuousBalancedClustering(G, 10);
    PartitionIntersection intersection;

    JaccardMeasure jaccard;

    EXPECT_EQ(1, jaccard.getDissimilarity(G, con, nonCon)) << "The Jaccard distance of a clustering with its complementary clustering should be 1";
    EXPECT_TRUE(GraphClusteringTools::isSingletonClustering(G, intersection.calculate(con, nonCon))) << "The intersection of a clustering with its complementary clustering should be the singleton clustering";
}

TEST_F(CommunityGTest, testHubDominance) {
    ClusteringGenerator generator;

    // make complete graph
    count n = 100;
    Graph G = Graph(n);
    G.forNodePairs([&](node u, node v) {
        G.addEdge(u,v);
    });
    Partition con = generator.makeContinuousBalancedClustering(G, 10);

    HubDominance hub;

    EXPECT_EQ(1u, hub.getQuality(con, G)) << "In complete graphs, the hub dominance is always 1";

    G = Graph(100);
    EXPECT_EQ(0u, hub.getQuality(con, G)) << "In graphs without nodes, the hub dominance is always 0";

    for (node v = 1; v < 4; ++v) {
        G.addEdge(0, v);
    }

    EXPECT_LT(0, hub.getQuality(con, G)) << "In graphs with internal edges, the hub dominance must be > 0";

    EXPECT_DOUBLE_EQ(3.0 / 9.0 * 1.0 / 10.0, hub.getQuality(con, G)) << "In this case, in one of ten equally-sized clusters a node has 3 of 9 possible neighbors.";

    Cover cov(con);

    EXPECT_DOUBLE_EQ(3.0 / 9.0 * 1.0 / 10.0, hub.getQuality(cov, G)) << "The Cover should have the same hub dominance as the equivalent partition";

    index newS = cov.upperBound();
    cov.setUpperBound(newS + 1);

    for (node u = 0; u < 10; ++u) {
        cov.addToSubset(newS, u);
    }

    EXPECT_DOUBLE_EQ(3.0 / 9.0 * 2.0 / 11.0, hub.getQuality(cov, G)) << "Duplicated subsets in covers count twice for the hub dominance";

    con = generator.makeSingletonClustering(G);

    EXPECT_EQ(1u, hub.getQuality(con, G)) << "The singleton partition has hub dominance 1 by definition.";

    cov = Cover(con);

    EXPECT_EQ(1u, hub.getQuality(con, G)) << "The singleton cover has hub dominance 1 by definition.";
}

TEST_F(CommunityGTest, testIntrapartitionDensity) {
    Aux::Random::setSeed(42, false);
    ClusteringGenerator generator;
    count n = 100;
    Graph G(n);
    G.forNodePairs([&](node u, node v) {
        G.addEdge(u, v);
    });
    Partition con(generator.makeContinuousBalancedClustering(G, 10));

    IntrapartitionDensity den(G, con);
    den.run();

    EXPECT_DOUBLE_EQ(1.0, den.getUnweightedAverage());
    EXPECT_DOUBLE_EQ(1.0, den.getMinimumValue());
    EXPECT_DOUBLE_EQ(1.0, den.getGlobal());

    G = Graph(100);

    den.run();

    EXPECT_DOUBLE_EQ(0.0, den.getUnweightedAverage());
    EXPECT_DOUBLE_EQ(0.0, den.getMinimumValue());
    EXPECT_DOUBLE_EQ(0.0, den.getGlobal());

    ClusteredRandomGraphGenerator gen(100, 2, 0.1, 0.02);
    G = gen.generate();
    auto P = gen.getCommunities();

    IntrapartitionDensity den2(G, P);
    den2.run();

    EXPECT_GT(1.0, den2.getUnweightedAverage());
    EXPECT_LT(0.0, den2.getUnweightedAverage());
    EXPECT_GT(1.0, den2.getMinimumValue());
    EXPECT_LT(0.0, den2.getMinimumValue());
    EXPECT_GT(1.0, den2.getGlobal());
    EXPECT_LT(0.0, den2.getGlobal());


    EXPECT_PRED_FORMAT2(::testing::DoubleLE, den.getMinimumValue(), den.getUnweightedAverage());
}

TEST_F(CommunityGTest, testPartitionFragmentation) {
    ClusteringGenerator generator;
    count n = 100;
    Graph G(n);
    G.forNodePairs([&](node u, node v) {
        G.addEdge(u, v);
    });
    Partition con(generator.makeContinuousBalancedClustering(G, 10));

    PartitionFragmentation frag1(G, con);
    frag1.run();
    EXPECT_EQ(0, frag1.getMaximumValue());
    EXPECT_EQ(0, frag1.getUnweightedAverage());
    EXPECT_EQ(0, frag1.getWeightedAverage());

    G = Graph(100);
    PartitionFragmentation frag2(G, con);
    frag2.run();
    EXPECT_DOUBLE_EQ(0.9, frag2.getMaximumValue());
    EXPECT_DOUBLE_EQ(0.9, frag2.getUnweightedAverage());
    EXPECT_DOUBLE_EQ(0.9, frag2.getWeightedAverage());

    con.setUpperBound(con.upperBound() + 1);
    G.forNodes([&](node u) {
        con[u] += 1;
    });

    frag2.run();

    EXPECT_DOUBLE_EQ(0.9, frag2.getMaximumValue());
    EXPECT_DOUBLE_EQ(0.9, frag2.getUnweightedAverage());
    EXPECT_DOUBLE_EQ(0.9, frag2.getWeightedAverage());

    PartitionFragmentation frag3(G, con);
    frag3.run();
    EXPECT_DOUBLE_EQ(0.9, frag3.getMaximumValue());
    EXPECT_DOUBLE_EQ(0.9, frag3.getUnweightedAverage());
    EXPECT_DOUBLE_EQ(0.9, frag3.getWeightedAverage());
}

TEST_F(CommunityGTest, testCoverF1Similarity) {
    count n = 20;
    Graph G(n);

    Cover C(n);
    C.setUpperBound(3);
    Cover ref(n);
    ref.setUpperBound(3);

    // 0: perfect overlap
    for (node u = 0; u < 10; ++u) {
        C.addToSubset(0, u);
        ref.addToSubset(1, u);
    }

    for (node u = 0; u < 11; ++u) {
        ref.addToSubset(0, u);
    }

    // 1: two node overlap with 0, one node overlap with 1
    for (node u = 9; u < 19; ++u) {
        C.addToSubset(1, u);
    }

    // 2: no overlap
    for (node u = 11; u < 20; ++u) {
        C.addToSubset(2, u);
    }

    CoverF1Similarity sim(G, C, ref);
    sim.run();

    EXPECT_DOUBLE_EQ(1.0, sim.getMaximumValue());
    EXPECT_DOUBLE_EQ(0.0, sim.getMinimumValue());
    EXPECT_DOUBLE_EQ(1.0, sim.getValue(0));
    const double pre = 2.0 / 11.0;
    const double re = 2.0 / 10.0;
    const double f1 = 2.0 * (pre * re) / (pre + re);
    EXPECT_DOUBLE_EQ(f1, sim.getValue(1));
    EXPECT_DOUBLE_EQ(0.0, sim.getValue(2));
    EXPECT_DOUBLE_EQ((1.0 + f1) / 3.0, sim.getUnweightedAverage());
    EXPECT_DOUBLE_EQ((1.0 * 10.0 + f1 * 10.0) / 29.0, sim.getWeightedAverage());
}


TEST_F(CommunityGTest, testLPPottsSequential) {
    count numClusters = 4;
    ClusteredRandomGraphGenerator gen(100, numClusters, 0.7, 0.03);
    Graph G = gen.generate();

    LPPotts algo(G);
    algo.run();
    Partition partition = algo.getPartition();

    EXPECT_LE(partition.numberOfSubsets(), numClusters + 1);
    EXPECT_GE(partition.numberOfSubsets(), numClusters);
}


TEST_F(CommunityGTest, testLPPottsParallel) {
    count numClusters = 4;
    ClusteredRandomGraphGenerator gen(100, numClusters, 0.7, 0.03);
    Graph G = gen.generate();

    LPPotts algo(G, 0.3, none, 20, true);
    algo.run();
    Partition partition = algo.getPartition();

    EXPECT_LE(partition.numberOfSubsets(), numClusters + 1);
    EXPECT_GE(partition.numberOfSubsets(), numClusters);
}


} /* namespace NetworKit */
