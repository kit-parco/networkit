/*
 * CoarseningGTest.cpp
 *
 *  Created on: 20.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <gtest/gtest.h>

#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/coarsening/ClusteringProjector.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/auxiliary/Log.hpp>

namespace NetworKit {

class CoarseningBenchmark: public testing::Test {};

TEST_F(CoarseningBenchmark, benchmarkCoarsening) {
    count n = 10000;
    count redF = 100; // reduction factor
    count k = n/redF;
    DEBUG("generating graph with ", n, " nodes");
    auto gen = ErdosRenyiGenerator(n, 0.05);
    Graph G = gen.generate();

    DEBUG("generating random partition");
    ClusteringGenerator clusteringGen;
    Partition zeta = clusteringGen.makeRandomClustering(G, k);

    //count k = zeta.numberOfSubsets();
    DEBUG("number of subsets: ", k);

    Aux::Timer timer;
    INFO("parallel coarsening");
    timer.start();
    ParallelPartitionCoarsening coarsening(G, zeta, true);
    coarsening.run();
    timer.stop();
    Graph Gc2 = coarsening.getCoarseGraph();
    INFO("parallel coarsening: ", timer.elapsedTag());
    EXPECT_EQ(k, Gc2.numberOfNodes());
}

} /* namespace NetworKit */

