// no-networkit-format
/*
 * GraphDistanceGTest.cpp
 *
 *  Created on: 23.07.2013
 *      Author: Henning Meyerhenke (meyerhenke@kit.edu)
 */

#include <gtest/gtest.h>

#include <networkit/distance/GraphDistance.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class GraphDistanceGTest: public testing::Test{};

// TODO: fix graph
TEST_F(GraphDistanceGTest, testGraphWeightedDistance) {
    METISGraphReader reader;
    Graph g = reader.read("input/grid-5x5-dist-arch.graph");

    GraphDistance gd;

    edgeweight dist_0_24 = gd.weightedDistance(g, 0, 24);
    edgeweight dist_24_0 = gd.weightedDistance(g, 24, 0);
    edgeweight dist_0_4 = gd.weightedDistance(g, 0, 4);
    edgeweight dist_4_0 = gd.weightedDistance(g, 4, 0);
    edgeweight dist_1_23 = gd.weightedDistance(g, 1, 23);
    edgeweight dist_5_19 = gd.weightedDistance(g, 5, 19);

    EXPECT_EQ(dist_0_24, 8);
    EXPECT_EQ(dist_24_0, 8);
    EXPECT_EQ(dist_0_4, 4);
    EXPECT_EQ(dist_4_0, 4);
    EXPECT_EQ(dist_1_23, 6);
    EXPECT_EQ(dist_5_19, 6);
}

TEST_F(GraphDistanceGTest, testGraphUnweightedDistance) {
    METISGraphReader reader;
    Graph g = reader.read("input/grid-5x5-dist-arch.graph");

    GraphDistance gd;

    edgeweight dist_0_24 = gd.unweightedDistance(g, 0, 24);
    edgeweight dist_24_0 = gd.unweightedDistance(g, 24, 0);
    edgeweight dist_0_4 = gd.unweightedDistance(g, 0, 4);
    edgeweight dist_4_0 = gd.unweightedDistance(g, 4, 0);
    edgeweight dist_1_23 = gd.unweightedDistance(g, 1, 23);
    edgeweight dist_5_19 = gd.unweightedDistance(g, 5, 19);

    EXPECT_EQ(dist_0_24, 8);
    EXPECT_EQ(dist_24_0, 8);
    EXPECT_EQ(dist_0_4, 4);
    EXPECT_EQ(dist_4_0, 4);
    EXPECT_EQ(dist_1_23, 6);
    EXPECT_EQ(dist_5_19, 6);
}


} /* namespace NetworKit */
