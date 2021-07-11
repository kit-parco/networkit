// no-networkit-format
/*
 * GlobalGTest.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include <gtest/gtest.h>

#include <networkit/global/ClusteringCoefficient.hpp>

#include <networkit/generators/ErdosRenyiGenerator.hpp>

namespace NetworKit {

class GlobalGTest: public testing::Test {};

TEST_F(GlobalGTest, testClusteringCoefficient) {

    ErdosRenyiGenerator graphGen(10, 1.0);
    Graph G = graphGen.generate();

    ClusteringCoefficient clusteringCoefficient;
    double cc = clusteringCoefficient.avgLocal(G);

    EXPECT_EQ(1.0, cc);
}


TEST_F(GlobalGTest, testGlobalClusteringCoefficient) {
    Graph G(6);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(1, 4);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(2, 5);
    G.addEdge(3, 5);

    double ccg = ClusteringCoefficient::exactGlobal(G);
    EXPECT_NEAR(ccg, 18.0 / 34.0, 1e-9);
}

} /* namespace NetworKit */
