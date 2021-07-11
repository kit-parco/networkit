// no-networkit-format
/*
 * DistMeasureTest.cpp
 *
 *  Created on: Jun 11, 2013
 *      Author: Henning
 */

#include <gtest/gtest.h>
#include <cstdio>

#include <networkit/graph/Graph.hpp>
#include <networkit/viz/PostscriptWriter.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/DibapGraphReader.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/linkprediction/AlgebraicDistanceIndex.hpp>

#include <networkit/auxiliary/Log.hpp>

namespace NetworKit {

class DistMeasuresGTest: public testing::Test {};

TEST_F(DistMeasuresGTest, testAlgebraicDistanceIndex) {
    Graph G(42);
    G.forNodePairs([&](node u, node v){
        G.addEdge(u,v);
    });

    count numSystems = 2;
    count numIterations = 200;
    double omega = 0.5;
    AlgebraicDistanceIndex ad(G, numSystems, numIterations, omega);
    ad.preprocess();

    double adSum = 0.0;
    G.forNodePairs([&](node u, node v){
        adSum += ad.run(u, v);
    });

    DEBUG("sum of algebraic distances: " , adSum);
    EXPECT_GE(1e-12, adSum) << "algebraic distances should converge towards zero";
}

} /* namespace NetworKit */

