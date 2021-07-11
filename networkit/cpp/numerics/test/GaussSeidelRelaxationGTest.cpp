// no-networkit-format
/*
 * GaussSeidelRelaxationGTest.cpp
 *
 *  Created on: 03.11.2014
 *      Author: Michael
 */
#include <gtest/gtest.h>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/numerics/GaussSeidelRelaxation.hpp>

namespace NetworKit {

class GaussSeidelRelaxationGTest : public testing::Test {};

TEST(GaussSeidelRelaxationGTest, debugSolve) {
    std::vector<Triplet> triplets = {{0,0,10}, {0,1,-1}, {0,2,2}, {1,0,-1}, {1,1,11}, {1,2,-1}, {1,3,3}, {2,0,2}, {2,1,-1}, {2,2,10}, {2,3,-1}, {3,1,3}, {3,2,-1}, {3,3,8}};
//	10  -1   2   0
//	-1  11  -1   3
//	 2  -1  10  -1
//	 0   3  -1   8
    CSRMatrix A(4, triplets);

    Vector b = {6, 25, -11, 15};
    Vector x = {0, 0, 0, 0};

    GaussSeidelRelaxation<CSRMatrix> solver;
    Vector result = solver.relax(A, b, x);

    EXPECT_EQ(1, std::round(result[0]));
    EXPECT_EQ(2, std::round(result[1]));
    EXPECT_EQ(-1, std::round(result[2]));
    EXPECT_EQ(1, std::round(result[3]));
}

TEST(GaussSeidelRelaxationGTest, debugIteration) {
    std::vector<Triplet> triplets = {{0,0,10}, {0,1,-1}, {0,2,2}, {1,0,-1}, {1,1,11}, {1,2,-1}, {1,3,3}, {2,0,2}, {2,1,-1}, {2,2,10}, {2,3,-1}, {3,1,3}, {3,2,-1}, {3,3,8}};
    //	10  -1   2   0
    //	-1  11  -1   3
    //	 2  -1  10  -1
    //	 0   3  -1   8
    CSRMatrix A(4, triplets);

    Vector b = {6, 25, -11, 15};
    Vector x = {0, 0, 0, 0};

    GaussSeidelRelaxation<CSRMatrix> solver;
    Vector result = solver.relax(A, b, x, 1);

    EXPECT_TRUE(result[0] > 0);
    EXPECT_TRUE(result[1] > 1);
    EXPECT_TRUE(result[2] < 0);
    EXPECT_TRUE(result[3] > 0);
}

} /* namespace NetworKit */
