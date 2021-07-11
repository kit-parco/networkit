/*
 * AllSimplePathsGTest.cpp
 *
 *  Created on: 27.06.2017
 *      Author: Eugenio Angriman
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <string>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/reachability/AllSimplePaths.hpp>

namespace NetworKit {

class AllSimplePathsGTest : public testing::Test {};

TEST_F(AllSimplePathsGTest, testAllSimplePaths) {
    EdgeListReader reader('\t', 0, "#", true, true);
    Graph G = reader.read("input/example.edgelist");

    AllSimplePaths allSimplePaths(G, 1, 9);
    EXPECT_ANY_THROW(allSimplePaths.run());

    G.addEdge(9, 6);
    G.addEdge(6, 9);
    AllSimplePaths allSimplePaths2(G, 3, 1);
    EXPECT_NO_THROW(allSimplePaths2.run());

    ASSERT_EQ(allSimplePaths2.numberOfSimplePaths(), 4);
    std::vector<node> path1{3, 7, 10, 9, 6, 1};
    std::vector<node> path2{3, 7, 10, 9, 6, 5, 1};
    std::vector<node> path3{3, 10, 9, 6, 1};
    std::vector<node> path4{3, 10, 9, 6, 5, 1};
    std::vector<std::vector<node>> results{path1, path2, path3, path4};

    allSimplePaths2.parallelForAllSimplePaths([&](const std::vector<node> &p) {
        ASSERT_TRUE(std::find(results.begin(), results.end(), p) != results.end());
    });

    std::vector<std::vector<node>> paths = allSimplePaths2.getAllSimplePaths();

    // backwards check
    for (count i = 0; i < results.size(); i++)
        ASSERT_TRUE(std::find(paths.begin(), paths.end(), results[i]) != paths.end());
}

} /* namespace NetworKit */
