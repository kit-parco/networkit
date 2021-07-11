// no-networkit-format
/*
 * OverlapGTest.cpp
 *
 *  Created on: 21.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <gtest/gtest.h>
#include <functional>

#include <networkit/overlap/HashingOverlapper.hpp>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/community/GraphClusteringTools.hpp>

namespace NetworKit {

class OverlapGTest: public testing::Test {};

TEST_F(OverlapGTest, testHashingOverlapperOnSingletonClusterings) {
    int64_t n = 10;
    Graph G(n);
    G.forNodePairs([&](node u, node v){
        G.addEdge(u,v);
    });

    ClusteringGenerator clusterGen;
    std::vector<Partition> clusterings;
    int z = 3; // number of clusterings
    for (int i = 0; i < z; ++i) {
        clusterings.push_back(clusterGen.makeSingletonClustering(G));
    }

    HashingOverlapper over;
    Partition core = over.run(G, clusterings);

    // test if core clustering is singleton-clustering
    bool isSingleton = true;
    G.forEdges([&](node u, node v) {
        isSingleton = isSingleton && (core.subsetOf(u) != core.subsetOf(v));
    });

    EXPECT_TRUE(isSingleton) << "overlap of multiple  singleton clusterings should be a singleton clustering";
}

TEST_F(OverlapGTest, testHashingOverlapperOnOneClusterings) {
    int64_t n = 10;
    Graph G(n);
    G.forNodePairs([&](node u, node v){
        G.addEdge(u,v);
    });

    ClusteringGenerator clusterGen;
    std::vector<Partition> clusterings;
    int z = 3; // number of clusterings
    for (int i = 0; i < z; ++i) {
        clusterings.push_back(clusterGen.makeOneClustering(G));
    }
    DEBUG("end of loop");

    HashingOverlapper over;
    Partition core = over.run(G, clusterings);

    // test if core clustering is one-clustering
    node v = 1;
    index one = core.subsetOf(v);
    bool isOneClustering = true; //TODO function call?
    G.forNodes([&](node v) {
        index c = core.subsetOf(v);
        DEBUG("CLUSTER! c = " , c);
        isOneClustering = isOneClustering && (c == one);
    });

    EXPECT_TRUE(isOneClustering) << "overlap of multiple 1-clustering should be a 1-clustering";

}

TEST_F(OverlapGTest, testHashingOverlapperForCorrectness) {
    count n = 4;
    Graph G(n);

    Partition zeta(n);
    Partition eta(n);

    zeta.setUpperBound(2);
    zeta[0] = 0;
    zeta[1] = 0;
    zeta[2] = 1;
    zeta[3] = 1;

    eta.setUpperBound(2);
    eta[0] = 0;
    eta[1] = 1;
    eta[2] = 0;
    eta[3] = 1;

    std::vector<Partition> clusterings = {zeta, eta};
    HashingOverlapper overlapper;
    Partition overlap = overlapper.run(G, clusterings);
    std::vector<node> overlapping_comparison = {0, 1, 2, 3};
    std::vector<node> overlapping_origin = overlap.getVector();
    EXPECT_EQ(overlap.numberOfSubsets(), 4);
    EXPECT_TRUE(std::is_permutation(overlapping_origin.begin(),overlapping_origin.end(),overlapping_comparison.begin()));
}

TEST_F(OverlapGTest, debugHashingOverlapperCorrectness) {
    count n = 2;
    Graph G(n);

    std::vector<Partition> clusterings;
    clusterings.emplace_back(n);
    clusterings[0].setUpperBound(10);
    clusterings[0][0] = 3;
    clusterings[0][1] = 9;
    clusterings.emplace_back(n);
    clusterings[1].setUpperBound(7);
    clusterings[1][0] = 6;
    clusterings[1][1] = 2;
    clusterings.emplace_back(n);
    clusterings[2].setUpperBound(1);
    clusterings[2][0] = 0;
    clusterings[2][1] = 0;

    HashingOverlapper overlapper;
    Partition overlap = overlapper.run(G, clusterings);

    EXPECT_TRUE(GraphClusteringTools::isSingletonClustering(G, overlap)) << "When one singleton clustering is in the overlap, the result should be a singleton clustering";
}

TEST_F(OverlapGTest, debugHashingOverlapperCorrectness2) {
    count n = 10000;
    count k = 1000;
    Graph G(n);
    ClusteringGenerator generator;

    std::vector<Partition> clusterings;
    clusterings.push_back(generator.makeRandomClustering(G, k));
    clusterings.push_back(generator.makeSingletonClustering(G));
    clusterings.push_back(generator.makeContinuousBalancedClustering(G, k));

    HashingOverlapper overlapper;
    Partition overlap = overlapper.run(G, clusterings);
    INFO("Number of clusters in the overlap is ", overlap.numberOfSubsets(), " and should be ", n);

    EXPECT_TRUE(GraphClusteringTools::isSingletonClustering(G, overlap)) << "When a singleton clustering is in the overlap, the result should be a singleton clustering";
}

} /* namespace NetworKit */
