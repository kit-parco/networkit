// no-networkit-format
/*
 * OctreeGTest.cpp
 *
 *  Created on: Apr 21, 2016
 *      Author: Michael
 */

#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/viz/Octree.hpp>

namespace NetworKit {

class OctreeGTest : public testing::Test {};

TEST(OctreeGTest, testOctreeWithExample) {
    std::vector<Vector> coordinates(2, Vector(11, 0.0));
    coordinates[0][0] = 2.0;
    coordinates[1][0] = 1.0;

    coordinates[0][1] = 1.0;
    coordinates[1][1] = 2.0;

    coordinates[0][2] = 7.0;
    coordinates[1][2] = 7.0;

    coordinates[0][3] = 8.0;
    coordinates[1][3] = 9.0;

    coordinates[0][4] = 10.0;
    coordinates[1][4] = 7.0;

    coordinates[0][5] = 2.0;
    coordinates[1][5] = 22.0;

    coordinates[0][6] = 10.0;
    coordinates[1][6] = 14.0;

    coordinates[0][7] = 15.0;
    coordinates[1][7] = 9.0;

    coordinates[0][8] = 19.0;
    coordinates[1][8] = 5.0;

    coordinates[0][9] = 23.0;
    coordinates[1][9] = 1.0;

    coordinates[0][10] = 21.0;
    coordinates[1][10] = 22.0;


    Octree<double> ocTree(coordinates);
    DEBUG(ocTree.toString());
    std::vector<std::pair<count, Point<double>>> fiveApprox = ocTree.approximateDistance(Point<double>(2.0, 22.0), 0.5);

    EXPECT_EQ(fiveApprox[0].first, 2u);
    EXPECT_EQ(fiveApprox[1].first, 3u);
    EXPECT_EQ(fiveApprox[2].first, 3u);
    EXPECT_EQ(fiveApprox[3].first, 1u);
    EXPECT_EQ(fiveApprox[4].first, 1u);

    for (count i = 0; i < 50; ++i) {
        Point<double> queryPoint = {(double)Aux::Random::integer(0, 24), (double)Aux::Random::integer(0, 24)};
        std::vector<std::pair<count, Point<double>>> result = ocTree.approximateDistance(queryPoint, 0.5);
        count sum = 0;
        for (auto &point : result) {
            sum += point.first;
        }

        EXPECT_NEAR(sum, 10.5, 0.5);
    }

    std::vector<std::pair<count, Point<double>>> exactApprox = ocTree.approximateDistance(Point<double>(2.0, 22.0), 0.0);
    EXPECT_EQ(exactApprox.size(), 10u);
}
} /* namespace NetworKit */
