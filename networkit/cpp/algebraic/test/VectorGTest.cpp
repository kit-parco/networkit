// no-networkit-format
/*
 * VectorGTest.cpp
 *
 *  Created on: 13.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */
#include <gtest/gtest.h>

#include <networkit/algebraic/Vector.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <cmath>

namespace NetworKit {


class VectorGTest : public testing::Test {};

TEST(VectorGTest, testVectorConstruction) {
    Vector v1;
    ASSERT_EQ(0u, v1.getDimension());

    Vector v2(10, 1.0);
    ASSERT_EQ(10u, v2.getDimension());
    for (uint64_t i = 0; i < v2.getDimension(); i++) {
        ASSERT_EQ(1.0, v2[i]);
    }

    Vector v3 = {1.0, 2.0, 3.0, 4.0, 5.0};
    ASSERT_EQ(5u, v3.getDimension());

    for (uint64_t i = 0; i < v3.getDimension(); i++) {
        EXPECT_EQ(static_cast<double>((i + 1)), v3[i]);
    }
}

TEST(VectorGTest, testDimension) {
    Vector testVector = {1.0, 2.0, 3.0, 4.0, 5.0};
    ASSERT_EQ(5u, testVector.getDimension());
}

TEST(VectorGTest, testTransposition) {
    Vector v = {1, 2, 3, 4, 5};
    EXPECT_FALSE(v.isTransposed());

    Vector v2 = v;
    EXPECT_FALSE(v2.isTransposed());

    v.transpose(); // returns transposed v, v should not be transposed
    EXPECT_FALSE(v.isTransposed());

    v2 = v.transpose();
    EXPECT_TRUE(v2.isTransposed());
}

TEST(VectorGTest, testLength) {
    Vector v = {12, 3, 9, 28, 0, -1};
    double result = 0.0;
    v.forElements([&](const double &value) {
        result += value * value;
    });

    EXPECT_EQ(std::sqrt(result), v.length());
}

TEST(VectorGTest, testAccessVectorElement) {
    Vector testVector = {1.0, 2.0, 3.0, 4.0, 5.0};

    double third = testVector[2];
    EXPECT_EQ(3.0, third);

    testVector[2]++;
    third = testVector[2];
    EXPECT_EQ(4.0, third);

    EXPECT_EQ(5.0, testVector.at(4));
}

TEST(VectorGTest, testOuterProduct) {
    Vector v1 = {1.0, -1.0};
    Vector v2 = {1.0, 2.0,  3.0};

    DynamicMatrix result = Vector::outerProduct(v1, v2);
    ASSERT_EQ(v1.getDimension(), result.numberOfRows());
    ASSERT_EQ(v2.getDimension(), result.numberOfColumns());

    EXPECT_EQ(1.0, result(0,0));
    EXPECT_EQ(2.0, result(0,1));
    EXPECT_EQ(3.0, result(0,2));
    EXPECT_EQ(-1.0, result(1,0));
    EXPECT_EQ(-2.0, result(1,1));
    EXPECT_EQ(-3.0, result(1,2));
}

TEST(VectorGTest, testInnerProduct) {
    Vector v1 = {1.0, 0.0, -1.0, -5.0, 2.0};
    Vector v2 = {1.0, 2.0,  3.0,  0.0, 5.0};
    Vector v3 = {1.0, 2.0};

    double dotProduct = v1.transpose() * v2;
    EXPECT_EQ(8.0, dotProduct);

    dotProduct = Vector::innerProduct(v1, v2);
    EXPECT_EQ(8.0, dotProduct);
}

TEST(VectorGTest, testVectorComparisonOperators) {
    Vector v1 = {1.0, 2.0, 3.0, 4.0, 5.0};
    Vector v2 = {1.0, 2.0, 3.0, 4.0, 5.0};

    Vector v3(5, 1.0);

    EXPECT_TRUE(v1 == v2);
    EXPECT_FALSE(v1 != v2);

    EXPECT_FALSE(v1 == v3);
    EXPECT_TRUE(v1 != v3);
}

TEST(VectorGTest, testVectorScalarMultiplication) {
    Vector testVector = {1.0, 2.0, 3.0, 4.0, 5.0};

    Vector vectorScalar = testVector * 2.0;
    ASSERT_EQ(testVector.getDimension(), vectorScalar.getDimension());

    Vector scalarVector = 2.0 * testVector;
    ASSERT_EQ(testVector.getDimension(), scalarVector.getDimension());

    for (uint64_t i = 0; i < testVector.getDimension(); i++) {
        EXPECT_EQ((i + 1) * 2.0, vectorScalar[i]);
        EXPECT_EQ((i + 1) * 2.0, scalarVector[i]);
    }
}

TEST(VectorGTest, testVectorMatrixMultiplication) {
    Vector v = {1.0, 2.0, 3.0};
    Vector v2 = {1.0, 2.0};
    // 8 3 4
    // 3 5 9
    // 4 9 2
    std::vector<Triplet> triplets = {{0,0,8}, {0,1,3}, {0,2,4}, {1,0,3}, {1,1,5}, {1,2,9}, {2,0,4}, {2,1,9}, {2,2,2}};

    DynamicMatrix mat(3,3,triplets);

    Vector result = v.transpose() * mat;
    ASSERT_TRUE(result.isTransposed());

    EXPECT_EQ(26, result[0]);
    EXPECT_EQ(40, result[1]);
    EXPECT_EQ(28, result[2]);
}

TEST(VectorGTest, testVectorDivisionOperator) {
    Vector testVector = {1.0, 2.0, 3.0, 4.0, 5.0};

    Vector vectorScalar = testVector / (1.0 / 2.0);
    ASSERT_EQ(testVector.getDimension(), vectorScalar.getDimension());

    for (uint64_t i = 0; i < testVector.getDimension(); i++) {
        EXPECT_EQ((i + 1) * 2.0, vectorScalar[i]);
    }
}

TEST(VectorGTest, testVectorAddition) {
    Vector v1 = {1.0, 2.0, 3.0, 4.0, 5.0};
    Vector v2 = {1.0, 2.0, 3.0, 4.0, 5.0};
    Vector v4 = {1.0, 2.0};

    Vector v3 = v1 + v2;
    ASSERT_EQ(v1.getDimension(), v3.getDimension());
    for (uint64_t i = 0; i < v3.getDimension(); i++) {
        EXPECT_EQ(v1[i] + v2[i], v3[i]);
    }

    v1 += v2;
    for (uint64_t i = 0; i < v1.getDimension(); i++) {
        EXPECT_EQ(v3[i], v1[i]);
    }

    v3 = v1.transpose() + v2.transpose();
    EXPECT_TRUE(v3.isTransposed());
}

TEST(VectorGTest, testVectorSubtraction) {
    Vector v1 = {2.0, 4.0, 6.0, 8.0, 10.0};
    Vector v2 = {1.0, 2.0, 3.0, 4.0, 5.0};
    Vector v4 = {1.0, 2.0};

    Vector v3 = v1 - v2;
    ASSERT_EQ(v1.getDimension(), v3.getDimension());
    for (uint64_t i = 0; i < v3.getDimension(); i++) {
        EXPECT_EQ(v1[i] - v2[i], v3[i]);
    }

    v1 -= v2;
    for (uint64_t i = 0; i < v1.getDimension(); i++) {
        EXPECT_EQ(v3[i], v1[i]);
    }

    v1 = Vector{1.0, 2.0, 3.0, 4.0};
    v1 -= 1.0;
    EXPECT_EQ(v1[0], 0.0);
    EXPECT_EQ(v1[1], 1.0);
    EXPECT_EQ(v1[2], 2.0);
    EXPECT_EQ(v1[3], 3.0);
}

TEST(VectorGTest, testVectorIterators) {
    Vector v = {1.0, 2.0, 3.0, 4.0, 5.0};

    // constant iterators
    int idx = 0;
    auto constantTester = [&](const double &element) {
        EXPECT_EQ(v[idx], element);
        idx++;
    };

    v.forElements(constantTester);

    auto constantParallelTester = [&](const int &idx, const double &element) {
        EXPECT_EQ(v[idx], element);
    };
    v.parallelForElements(constantParallelTester);

    // non-constant iterators
    auto nonConstantTester = [&](double &element) {
        element++;
    };

    v.forElements(nonConstantTester);
    for (uint64_t i = 0; i < v.getDimension(); i++) {
        EXPECT_EQ((i + 2.0), v[i]);
    }

    auto nonConstantParallelTester = [](const int &, double &element) {
        element++;
    };

    v.parallelForElements(nonConstantParallelTester);
    for (uint64_t i = 0; i < v.getDimension(); i++) {
        EXPECT_EQ((i + 3.0), v[i]);
    }
}

TEST(VectorGTest, testMean) {
    Vector v = {1.0, 2.0, 3.0, 4.0, 5.0};
    double mean = v.mean();
    EXPECT_EQ(mean, 3);
}


} /* namespace NetworKit */
