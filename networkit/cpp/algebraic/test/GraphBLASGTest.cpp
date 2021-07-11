// no-networkit-format
/*
 * GraphBLASGTest.cpp
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <gtest/gtest.h>
#include <iostream>

#include <networkit/algebraic/GraphBLAS.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>

namespace NetworKit {

class GraphBLASGTest : public testing::Test {};

TEST_F(GraphBLASGTest, testMxM) {
    std::vector<Triplet> triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,2,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};

    //
    //				 1  2  3  0
    // 				 2  2  0  0
    // mat1 = mat2 = 3  0  3 -1
    //				 0  0 -1  4
    //
    CSRMatrix mat1(4, triplets);
    CSRMatrix mat2(4, triplets);

    //
    //			14  6  12  -3
    //			 6  8   6   0
    // result = 12  6  19  -7
    //			-3  0  -7  17
    //
    CSRMatrix result = GraphBLAS::MxM(mat1, mat2);

    EXPECT_EQ(14, result(0,0));
    EXPECT_EQ(6, result(0,1));
    EXPECT_EQ(12, result(0,2));
    EXPECT_EQ(-3, result(0,3));
    EXPECT_EQ(6, result(1,0));
    EXPECT_EQ(8, result(1,1));
    EXPECT_EQ(6, result(1,2));
    EXPECT_EQ(0, result(1,3));
    EXPECT_EQ(12, result(2,0));
    EXPECT_EQ(6, result(2,1));
    EXPECT_EQ(19, result(2,2));
    EXPECT_EQ(-7, result(2,3));
    EXPECT_EQ(-3, result(3,0));
    EXPECT_EQ(0, result(3,1));
    EXPECT_EQ(-7, result(3,2));
    EXPECT_EQ(17, result(3,3));


    // max-plus semiring
    mat1 = CSRMatrix(4, triplets, MaxPlusSemiring::zero());
    mat2 = CSRMatrix(4, triplets, MaxPlusSemiring::zero());

    //
    //			6  4  6  2
    //			4  4  5  zero()
    // result = 6  5  6  3
    //			2  zero()  3  8
    //
    result = GraphBLAS::MxM<MaxPlusSemiring>(mat1, mat2);

    EXPECT_EQ(6, result(0,0));
    EXPECT_EQ(4, result(0,1));
    EXPECT_EQ(6, result(0,2));
    EXPECT_EQ(2, result(0,3));
    EXPECT_EQ(4, result(1,0));
    EXPECT_EQ(4, result(1,1));
    EXPECT_EQ(5, result(1,2));
    EXPECT_EQ(MaxPlusSemiring::zero(), result(1,3));
    EXPECT_EQ(6, result(2,0));
    EXPECT_EQ(5, result(2,1));
    EXPECT_EQ(6, result(2,2));
    EXPECT_EQ(3, result(2,3));
    EXPECT_EQ(2, result(3,0));
    EXPECT_EQ(MaxPlusSemiring::zero(), result(3,1));
    EXPECT_EQ(3, result(3,2));
    EXPECT_EQ(8, result(3,3));

    triplets[7] = {2,3,1};
    triplets[8] = {3,2,1};
    //
    //				 1  2  3  0
    // 				 2  2  0  0
    // mat1 = mat2 = 3  0  3  1
    //				 0  0  1  4
    //
    mat1 = CSRMatrix(4, triplets, MinMaxSemiring::zero());
    mat2 = CSRMatrix(4, triplets, MinMaxSemiring::zero());


    result = GraphBLAS::MxM<MinMaxSemiring>(mat1, mat2);

    EXPECT_EQ(1, result(0,0));
    EXPECT_EQ(2, result(0,1));
    EXPECT_EQ(3, result(0,2));
    EXPECT_EQ(3, result(0,3));
    EXPECT_EQ(2, result(1,0));
    EXPECT_EQ(2, result(1,1));
    EXPECT_EQ(3, result(1,2));
    EXPECT_EQ(MinMaxSemiring::zero(), result(1,3));
    EXPECT_EQ(3, result(2,0));
    EXPECT_EQ(3, result(2,1));
    EXPECT_EQ(1, result(2,2));
    EXPECT_EQ(3, result(2,3));
    EXPECT_EQ(3, result(3,0));
    EXPECT_EQ(MinMaxSemiring::zero(), result(3,1));
    EXPECT_EQ(3, result(3,2));
    EXPECT_EQ(1, result(3,3));
}

TEST_F(GraphBLASGTest, testMxMAccum) {
    std::vector<Triplet> triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,2,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};

    //
    //				 1  2  3  0
    // 				 2  2  0  0
    // mat1 = mat2 = 3  0  3 -1
    //				 0  0 -1  4
    //
    CSRMatrix mat1(4, triplets);
    CSRMatrix mat2(4, triplets);

    CSRMatrix C = CSRMatrix::diagonalMatrix(Vector({-2,1,3,0}));

    //
    //			12  6  12  -3
    //			 6  9   6   0
    // result = 12  6  22  -7
    //			-3  0  -7  17
    //
    GraphBLAS::MxM(mat1, mat2, C);

    EXPECT_EQ(12, C(0,0));
    EXPECT_EQ(6, C(0,1));
    EXPECT_EQ(12, C(0,2));
    EXPECT_EQ(-3, C(0,3));
    EXPECT_EQ(6, C(1,0));
    EXPECT_EQ(9, C(1,1));
    EXPECT_EQ(6, C(1,2));
    EXPECT_EQ(0, C(1,3));
    EXPECT_EQ(12, C(2,0));
    EXPECT_EQ(6, C(2,1));
    EXPECT_EQ(22, C(2,2));
    EXPECT_EQ(-7, C(2,3));
    EXPECT_EQ(-3, C(3,0));
    EXPECT_EQ(0, C(3,1));
    EXPECT_EQ(-7, C(3,2));
    EXPECT_EQ(17, C(3,3));


    GraphBLAS::MxM(mat1, mat2, C, [](const double, const double b){return b;});

    EXPECT_EQ(14, C(0,0));
    EXPECT_EQ(6, C(0,1));
    EXPECT_EQ(12, C(0,2));
    EXPECT_EQ(-3, C(0,3));
    EXPECT_EQ(6, C(1,0));
    EXPECT_EQ(8, C(1,1));
    EXPECT_EQ(6, C(1,2));
    EXPECT_EQ(0, C(1,3));
    EXPECT_EQ(12, C(2,0));
    EXPECT_EQ(6, C(2,1));
    EXPECT_EQ(19, C(2,2));
    EXPECT_EQ(-7, C(2,3));
    EXPECT_EQ(-3, C(3,0));
    EXPECT_EQ(0, C(3,1));
    EXPECT_EQ(-7, C(3,2));
    EXPECT_EQ(17, C(3,3));
}

TEST_F(GraphBLASGTest, testMxV) {
    std::vector<Triplet> triplets;

    for (index i = 0; i < 10000; ++i) {
        triplets.push_back({i,i, (double) i});
    }

    triplets.push_back({42,43,42.0});

    Vector vector(10000, 1.0);
    vector[500] = 3.5;

    CSRMatrix mat(10000, triplets);

    Vector result = GraphBLAS::MxV(mat, vector);
    ASSERT_EQ(mat.numberOfRows(), result.getDimension());

    for (index i = 0; i < 1000; ++i) {
        if (i != 500 && i != 42 && i != 43) {
            EXPECT_EQ(i, result[i]);
        }
    }

    EXPECT_EQ(42.0, mat(42, 43));
    EXPECT_EQ(84.0, result[42]);
    EXPECT_EQ(1750.0, result[500]);

    //		  1  2  3  0
    //        2  2  0  0
    // mat2 = 3  0  3 -1
    //		  0  0 -1  4
    triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,2,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};

    CSRMatrix mat2(4, triplets);

    Vector v({1,2,3,0});
    Vector res = GraphBLAS::MxV(mat2, v);
    ASSERT_EQ(mat2.numberOfRows(), res.getDimension());

    EXPECT_EQ(14, res[0]);
    EXPECT_EQ(6, res[1]);
    EXPECT_EQ(12, res[2]);
    EXPECT_EQ(-3, res[3]);

    v = {1,2,3,MinPlusSemiring::zero()};
    mat2 = CSRMatrix(4, triplets, MinPlusSemiring::zero());
    res = GraphBLAS::MxV<MinPlusSemiring>(mat2, v);

    EXPECT_EQ(2, res[0]);
    EXPECT_EQ(3, res[1]);
    EXPECT_EQ(4, res[2]);
    EXPECT_EQ(2, res[3]);
}

TEST_F(GraphBLASGTest, testEWiseAdd) {
    std::vector<Triplet> triplets1;
    std::vector<Triplet> triplets2;

    for (index i = 0; i < 1000; ++i) {
        triplets1.push_back({i,i,1});
        triplets2.push_back({i,i,(double)i});
    }

    triplets1.push_back({2,71, 1.8});
    triplets2.push_back({42,43,3.14});

    CSRMatrix mat1(1000, triplets1);
    CSRMatrix mat2(1000, triplets2);

    CSRMatrix result = GraphBLAS::eWiseAdd(mat1, mat2);
    ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

    EXPECT_EQ(0.0, result(10, 13));


    for (index i = 0; i < result.numberOfRows(); ++i) {
        EXPECT_EQ((i + 1), result(i, i));
    }
    EXPECT_EQ(1.8, result(2, 71));
    EXPECT_EQ(3.14, result(42, 43));
    EXPECT_EQ(0.0, result(3, 14));

    //		  1  2  3  0
    //        2  2  0  0
    // mat1 = 3  0  3 -1
    //		  0  0 -1  4
    std::vector<Triplet> triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,2,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};

    mat1 = CSRMatrix(4, triplets, MinMaxSemiring::zero());

    triplets[7] = {2,3,1};
    triplets[8] = {3,2,1};
    //
    //		  1  2  3  0
    // 	 	  2  2  0  0
    // mat2 = 3  0  3  1
    //		  0  0  1  4
    //
    mat2 = CSRMatrix(4, triplets, MinMaxSemiring::zero());

    result = GraphBLAS::eWiseAdd<MinMaxSemiring>(mat1, mat2);
    EXPECT_EQ(1, result(0,0));
    EXPECT_EQ(2, result(0,1));
    EXPECT_EQ(3, result(0,2));
    EXPECT_EQ(MinMaxSemiring::zero(), result(0,3));
    EXPECT_EQ(2, result(1,0));
    EXPECT_EQ(2, result(1,1));
    EXPECT_EQ(MinMaxSemiring::zero(), result(1,2));
    EXPECT_EQ(MinMaxSemiring::zero(), result(1,3));
    EXPECT_EQ(3, result(2,0));
    EXPECT_EQ(MinMaxSemiring::zero(), result(2,1));
    EXPECT_EQ(3, result(2,2));
    EXPECT_EQ(-1, result(2,3));
    EXPECT_EQ(MinMaxSemiring::zero(), result(3,0));
    EXPECT_EQ(MinMaxSemiring::zero(), result(3,1));
    EXPECT_EQ(-1, result(3,2));
    EXPECT_EQ(4, result(3,3));
}

TEST_F(GraphBLASGTest, testEWiseMult) {
    std::vector<Triplet> triplets1;
    std::vector<Triplet> triplets2;

    for (index i = 0; i < 1000; ++i) {
        triplets1.push_back({i,i,1});
        triplets2.push_back({i,i,(double)i});
    }

    triplets1.push_back({2,71, 1.8});
    triplets2.push_back({42,43,3.14});

    CSRMatrix mat1(1000, triplets1);
    CSRMatrix mat2(1000, triplets2);

    CSRMatrix result = GraphBLAS::eWiseMult(mat1, mat2);
    ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

    EXPECT_EQ(0.0, result(10, 13));

    for (index i = 0; i < result.numberOfRows(); ++i) {
        EXPECT_EQ(i, result(i, i));
    }

    EXPECT_EQ(0.0, result(2, 71));
    EXPECT_EQ(0.0, result(42, 43));
    EXPECT_EQ(0.0, result(3, 14));


    //		  1  2  3  0
    //        2  2  0  0
    // mat1 = 3  0  3 -1
    //		  0  0 -1  4
    std::vector<Triplet> triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,2,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};

    mat1 = CSRMatrix(4, triplets);
    mat2 = CSRMatrix(4, triplets);

    result = GraphBLAS::eWiseMult(mat1, mat2);
    EXPECT_EQ(1, result(0,0));
    EXPECT_EQ(4, result(0,1));
    EXPECT_EQ(9, result(0,2));
    EXPECT_EQ(0, result(0,3));
    EXPECT_EQ(4, result(1,0));
    EXPECT_EQ(4, result(1,1));
    EXPECT_EQ(0, result(1,2));
    EXPECT_EQ(0, result(1,3));
    EXPECT_EQ(9, result(2,0));
    EXPECT_EQ(0, result(2,1));
    EXPECT_EQ(9, result(2,2));
    EXPECT_EQ(1, result(2,3));
    EXPECT_EQ(0, result(3,0));
    EXPECT_EQ(0, result(3,1));
    EXPECT_EQ(1, result(3,2));
    EXPECT_EQ(16, result(3,3));

    mat1 = CSRMatrix(4, triplets, MinMaxSemiring::zero());

    triplets[7] = {2,3,1};
    triplets[8] = {3,2,1};
    //
    //		  1  2  3  0
    // 	 	  2  2  0  0
    // mat2 = 3  0  3  1
    //		  0  0  1  4
    //
    mat2 = CSRMatrix(4, triplets, MinMaxSemiring::zero());

    result = GraphBLAS::eWiseMult<MinMaxSemiring>(mat1, mat2);
    EXPECT_EQ(1, result(0,0));
    EXPECT_EQ(2, result(0,1));
    EXPECT_EQ(3, result(0,2));
    EXPECT_EQ(MinMaxSemiring::zero(), result(0,3));
    EXPECT_EQ(2, result(1,0));
    EXPECT_EQ(2, result(1,1));
    EXPECT_EQ(MinMaxSemiring::zero(), result(1,2));
    EXPECT_EQ(MinMaxSemiring::zero(), result(1,3));
    EXPECT_EQ(3, result(2,0));
    EXPECT_EQ(MinMaxSemiring::zero(), result(2,1));
    EXPECT_EQ(3, result(2,2));
    EXPECT_EQ(1, result(2,3));
    EXPECT_EQ(MinMaxSemiring::zero(), result(3,0));
    EXPECT_EQ(MinMaxSemiring::zero(), result(3,1));
    EXPECT_EQ(1, result(3,2));
    EXPECT_EQ(4, result(3,3));
}

TEST_F(GraphBLASGTest, rowReductionTest) {
    //		  1  2  3  0
    //        2  2  0  0
    // mat1 = 3  0  3 -1
    //		  0  0 -1  4
    std::vector<Triplet> triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,2,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};
    CSRMatrix mat(4, triplets);

    Vector rowReduction = GraphBLAS::rowReduce(mat);
    ASSERT_EQ(rowReduction.getDimension(), mat.numberOfRows());
    EXPECT_EQ(6, rowReduction[0]);
    EXPECT_EQ(4, rowReduction[1]);
    EXPECT_EQ(5, rowReduction[2]);
    EXPECT_EQ(3, rowReduction[3]);

    mat = CSRMatrix(4, triplets, MaxPlusSemiring::zero());
    rowReduction = GraphBLAS::columnReduce<MaxPlusSemiring>(mat);
    ASSERT_EQ(rowReduction.getDimension(), mat.numberOfRows());
    EXPECT_EQ(3, rowReduction[0]);
    EXPECT_EQ(2, rowReduction[1]);
    EXPECT_EQ(3, rowReduction[2]);
    EXPECT_EQ(4, rowReduction[3]);
}

TEST_F(GraphBLASGTest, columnReductionTest) {
    //		  1  2  3  0
    //        2  2  0  0
    // mat1 = 3  0  3 -1
    //		  0  0 -1  4
    std::vector<Triplet> triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,2,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};
    CSRMatrix mat(4, triplets);

    Vector columnReduction = GraphBLAS::columnReduce(mat);
    ASSERT_EQ(columnReduction.getDimension(), mat.numberOfColumns());
    EXPECT_EQ(6, columnReduction[0]);
    EXPECT_EQ(4, columnReduction[1]);
    EXPECT_EQ(5, columnReduction[2]);
    EXPECT_EQ(3, columnReduction[3]);

    mat = CSRMatrix(4, triplets, MinMaxSemiring::zero());
    columnReduction = GraphBLAS::columnReduce<MinMaxSemiring>(mat);
    ASSERT_EQ(columnReduction.getDimension(), mat.numberOfColumns());
    EXPECT_EQ(0, columnReduction[0]);
    EXPECT_EQ(0, columnReduction[1]);
    EXPECT_EQ(-1, columnReduction[2]);
    EXPECT_EQ(-1, columnReduction[3]);
}

} /* namespace NetworKit */
