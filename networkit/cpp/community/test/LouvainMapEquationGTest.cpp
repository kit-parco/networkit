/*
 * LouvainMapEquationGTest.cpp
 *
 * Created on: 2019-10-30
 * Author: Armin Wiebigke
  */

#include <gtest/gtest.h>

#include <networkit/community/LouvainMapEquation.hpp>
#include <networkit/generators/ClusteredRandomGraphGenerator.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class MapEquationGTest : public testing::Test {
public:
    void SetUp() {
        Aux::Random::setSeed(435913, false);
    }
};

void addClique(Graph &graph, Partition &groundTruth, node lowestId, node highestId) {
    index subsetId = groundTruth.upperBound();
    groundTruth.setUpperBound(subsetId + 1);
    for (node i = lowestId; i <= highestId; ++i) {
        groundTruth.addToSubset(subsetId, i);
        for (node j = i + 1; j <= highestId; ++j) {
            graph.addEdge(i, j);
        }
    }
}

TEST_F(MapEquationGTest, testLocalMoveSmall) {
    Aux::Random::setSeed(2342556, false);
    Graph G(10);
    Partition groundTruth(10);
    addClique(G, groundTruth, 0, 4);
    addClique(G, groundTruth, 5, 9);
    G.addEdge(0, 9);
    G.addEdge(1, 8);
    G.addEdge(2, 7);

    LouvainMapEquation mapequation(G, false);
    mapequation.run();
    auto partition = mapequation.getPartition();

    EXPECT_EQ(partition.getSubsets(), groundTruth.getSubsets());
}

TEST_F(MapEquationGTest, testLocalMove) {
    Aux::Random::setSeed(2342556, true);
    ClusteredRandomGraphGenerator generator(100, 4, 0.5, 0.05);
    Graph G = generator.generate();
    Partition groundTruth = generator.getCommunities();

    LouvainMapEquation algo(G, false);
    algo.run();
    auto partition = algo.getPartition();

    EXPECT_EQ(partition.getSubsets(), groundTruth.getSubsets());
}

TEST_F(MapEquationGTest, testLocalMoveLargeHierarchical) {
    Aux::Random::setSeed(2342556, false);
    ClusteredRandomGraphGenerator generator(100, 4, 0.5, 0.05);
    Graph G = generator.generate();
    Partition groundTruth = generator.getCommunities();

    LouvainMapEquation algo(G, true);
    algo.run();
    auto partition = algo.getPartition();

    EXPECT_EQ(partition.getSubsets(), groundTruth.getSubsets());
}

}
