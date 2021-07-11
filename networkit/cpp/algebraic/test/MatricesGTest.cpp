// no-networkit-format
/*
 * MatricesGTest.cpp
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/algebraic/MatrixTools.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>

namespace NetworKit {

class MatricesGTest : public testing::Test {
protected:
    virtual void SetUp() {
        METISGraphReader reader;
        graph = reader.read("input/PGPgiantcompo.graph");
    }

    Graph graph;

public:
    MatricesGTest() = default;
    virtual ~MatricesGTest() = default;

    template<class Matrix>
    void testDimension();

    template<class Matrix>
    void testNNZInRow();

    template<class Matrix>
    void testRowAndColumnAccess();

    template<class Matrix>
    void testDiagonalVector();

    template<class Matrix>
    void testTranspose();

    template<class Matrix>
    void testExtract();

    template<class Matrix>
    void testAssign();

    template<class Matrix>
    void testApply();

    template<class Matrix>
    void testMatrixAddition();

    template<class Matrix>
    void testMatrixSubtraction();

    template<class Matrix>
    void testScalarMultiplication();

    template<class Matrix>
    void testMatrixDivisionOperator();

    template<class Matrix>
    void testMatrixVectorProduct();

    template<class Matrix>
    void testMatrixMultiplication();

    template<class Matrix>
    void testBigMatrixMultiplication();

    template<class Matrix>
    void testAdjacencyMatrix();

    template<class Matrix>
    void testDiagonalMatrix();

    template<class Matrix>
    void testIncidenceMatrix();

    template<class Matrix>
    void testLaplacianOfGraph();

    // TODO: Test other matrix classes

    // TODO: Test mmT multiplication, etc.!
};

template<class Matrix>
void MatricesGTest::testDimension() {
    Matrix mat(10);

    ASSERT_EQ(10u, mat.numberOfRows());
    ASSERT_EQ(10u, mat.numberOfColumns());

    mat = Matrix(5, 10, 0.0);
    ASSERT_EQ(5u, mat.numberOfRows());
    ASSERT_EQ(10u, mat.numberOfColumns());

    mat = Matrix(10, 5, 0.0);
    ASSERT_EQ(10u, mat.numberOfRows());
    ASSERT_EQ(5u, mat.numberOfColumns());
}

template<class Matrix>
void MatricesGTest::testNNZInRow() {
    /*
     * 1.0  0.0  2.0
     * 4.0  0.0  0.0
     * 0.0  0.0  0.0
     * 0.0  0.0  2.0
     */
    std::vector<Triplet> triplets = {{0,0,1.0}, {0,2,2.0}, {1,0,4.0}, {3,3,2.0}};

    Matrix mat(4, triplets);
    EXPECT_EQ(2u, mat.nnzInRow(0));
    EXPECT_EQ(1u, mat.nnzInRow(1));
    EXPECT_EQ(0u, mat.nnzInRow(2));
    EXPECT_EQ(1u, mat.nnzInRow(3));
}

template<class Matrix>
void MatricesGTest::testRowAndColumnAccess() {
    std::vector<Triplet> triplets;

    for (index i = 0; i < 1000; ++i) {
        triplets.push_back({3,i,(double)i});
    }

    triplets.push_back({10,10,42.123});

    Matrix mat(1000, triplets);

    Vector v = mat.row(3);
    ASSERT_EQ(mat.numberOfColumns(), v.getDimension());

    for (index i = 0; i < 1000; ++i) {
        EXPECT_EQ(i, v[i]);
    }

    v = mat.row(10);
    ASSERT_EQ(v.getDimension(), mat.numberOfColumns());
    ASSERT_TRUE(v.isTransposed());
    EXPECT_EQ(42.123, v[10]);

    v = mat.column(10);
    ASSERT_EQ(mat.numberOfRows(), v.getDimension());
    ASSERT_FALSE(v.isTransposed());

    EXPECT_EQ(10.0, v[3]);
    EXPECT_EQ(42.123, v[10]);


    // rectangular matrix
    // n x m (n < m)
    triplets.clear();
    triplets.push_back({4,9,11});
    triplets.push_back({0,0,42});

    mat = Matrix(5, 10, triplets);
    v = mat.row(0);
    ASSERT_EQ(v.getDimension(), 10u);
    for (index i = 0; i < v.getDimension(); ++i) {
        if (i == 0) {
            EXPECT_EQ(42, v[i]);
        } else {
            EXPECT_EQ(0, v[i]);
        }
    }

    v = mat.column(9);
    ASSERT_EQ(v.getDimension(), 5u);
    for (index i = 0; i < v.getDimension(); ++i) {
        if (i == v.getDimension() - 1) {
            EXPECT_EQ(11, v[i]);
        } else {
            EXPECT_EQ(0, v[i]);
        }
    }

    // rectangular matrix
    // n x m (n > m)

    triplets.clear();
    triplets.push_back({9,4,11});
    triplets.push_back({0,0,42});

    mat = Matrix(10, 5, triplets);
    v = mat.row(0);
    ASSERT_EQ(v.getDimension(), 5u);
    for (index i = 0; i < v.getDimension(); ++i) {
        if (i == 0) {
            EXPECT_EQ(42, v[i]);
        } else {
            EXPECT_EQ(0, v[i]);
        }
    }

    v = mat.column(4);
    ASSERT_EQ(v.getDimension(), 10u);
    for (index i = 0; i < v.getDimension(); ++i) {
        if (i == v.getDimension() - 1) {
            EXPECT_EQ(11, v[i]);
        } else {
            EXPECT_EQ(0, v[i]);
        }
    }
}

template<class Matrix>
void MatricesGTest::testDiagonalVector() {
    // 1  2  3  0
    // 2  2  0  0
    // 3  0  0 -1
    // 0  0 -1  4
    std::vector<Triplet> triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};

    Matrix mat(4,4,triplets);
    Vector diagonal = mat.diagonal();
    EXPECT_EQ(4u, diagonal.getDimension());
    EXPECT_EQ(1, diagonal[0]);
    EXPECT_EQ(2, diagonal[1]);
    EXPECT_EQ(0, diagonal[2]);
    EXPECT_EQ(4, diagonal[3]);

    // rectangular matrix
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    triplets = {{0,0,1}, {1,2,3}};
    mat = Matrix(2,5,triplets);
    diagonal = mat.diagonal();
    EXPECT_EQ(2u, diagonal.getDimension());
    EXPECT_EQ(1, diagonal[0]);
    EXPECT_EQ(0, diagonal[1]);

    // rectangular matrix
    //
    // 1 0
    // 0 -3
    // 0 0
    // 0 0
    // 0 0
    triplets = {{0,0,1.0}, {1,1,-3.0}};
    mat = Matrix(5,2,triplets);
    diagonal = mat.diagonal();
    EXPECT_EQ(2u, diagonal.getDimension());
    EXPECT_EQ(1, diagonal[0]);
    EXPECT_EQ(-3, diagonal[1]);
}

template<class Matrix>
void MatricesGTest::testTranspose() {
    // 1  0  1  0
    // 2  2  0  0
    // 3  0  0 -1
    // 0  0  1  4
    std::vector<Triplet> triplets = {{0,0,1}, {0,2,1}, {1,0,2}, {1,1,2}, {2,0,3}, {2,3,-1}, {3,2,1}, {3,3,4}};

    Matrix mat(4,4,triplets);
    Matrix matT = mat.transpose();
    EXPECT_EQ(4u, matT.numberOfRows());
    EXPECT_EQ(4u, matT.numberOfColumns());

    EXPECT_EQ(1, matT(0,0));
    EXPECT_EQ(2, matT(0,1));
    EXPECT_EQ(3, matT(0,2));
    EXPECT_EQ(0, matT(0,3));
    EXPECT_EQ(0, matT(1,0));
    EXPECT_EQ(2, matT(1,1));
    EXPECT_EQ(0, matT(1,2));
    EXPECT_EQ(0, matT(1,3));
    EXPECT_EQ(1, matT(2,0));
    EXPECT_EQ(0, matT(2,1));
    EXPECT_EQ(0, matT(2,2));
    EXPECT_EQ(1, matT(2,3));
    EXPECT_EQ(0, matT(3,0));
    EXPECT_EQ(0, matT(3,1));
    EXPECT_EQ(-1, matT(3,2));
    EXPECT_EQ(4, matT(3,3));

    // rectangular matrix
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    triplets = {{0,0,1}, {1,2,3}};
    mat = Matrix(2,5,triplets);
    matT = mat.transpose();
    EXPECT_EQ(5u, matT.numberOfRows());
    EXPECT_EQ(2u, matT.numberOfColumns());

    EXPECT_EQ(1, matT(0,0));
    EXPECT_EQ(0, matT(0,1));
    EXPECT_EQ(0, matT(1,0));
    EXPECT_EQ(0, matT(1,1));
    EXPECT_EQ(0, matT(2,0));
    EXPECT_EQ(3, matT(2,1));
    EXPECT_EQ(0, matT(3,0));
    EXPECT_EQ(0, matT(3,1));
    EXPECT_EQ(0, matT(4,0));
    EXPECT_EQ(0, matT(4,1));

    // rectangular matrix
    //
    // 1 0
    // 0 0
    // 0 3
    // 0 0
    // 0 0
    triplets = {{0,0,1.0}, {2,1,3.0}};
    mat = Matrix(5,2,triplets);
    matT = mat.transpose();
    EXPECT_EQ(2u, matT.numberOfRows());
    EXPECT_EQ(5u, matT.numberOfColumns());

    EXPECT_EQ(1, matT(0,0));
    EXPECT_EQ(0, matT(0,1));
    EXPECT_EQ(0, matT(0,2));
    EXPECT_EQ(0, matT(0,3));
    EXPECT_EQ(0, matT(0,4));
    EXPECT_EQ(0, matT(1,0));
    EXPECT_EQ(0, matT(1,1));
    EXPECT_EQ(3, matT(1,2));
    EXPECT_EQ(0, matT(1,3));
    EXPECT_EQ(0, matT(1,4));
}

template<class Matrix>
void MatricesGTest::testExtract() {
    Matrix mat = Matrix::adjacencyMatrix(graph);
    std::vector<index> rows(500);
    std::vector<index> columns(500);

    for (index i = 0; i < 500; ++i) {
        rows[i] = Aux::Random::integer(graph.numberOfNodes()-1);
        columns[i] = Aux::Random::integer(graph.numberOfNodes()-1);
    }

    Matrix subMatrix = mat.extract(rows, columns);
    ASSERT_EQ(rows.size(), subMatrix.numberOfRows());
    ASSERT_EQ(columns.size(), subMatrix.numberOfColumns());

    for (index i = 0; i < subMatrix.numberOfRows(); ++i) {
        for (index j = 0; j < subMatrix.numberOfColumns(); ++j) {
            EXPECT_EQ(mat(rows[i], columns[j]), subMatrix(i,j));
        }
    }


    // 1  0  1  0
    // 2  2  0  0
    // 3  0  0 -1
    // 0  0  1  4
    std::vector<Triplet> triplets = {{0,0,1}, {0,2,1}, {1,0,2}, {1,1,2}, {2,0,3}, {2,3,-1}, {3,2,1}, {3,3,4}};

    mat = Matrix(4,4,triplets);
    rows = {0,2,0};
    columns = {1,0,2,1};

    //
    // 0 1 1 0
    // 0 3 0 0
    // 0 1 1 0
    //
    subMatrix = mat.extract(rows, columns);
    ASSERT_EQ(3u, subMatrix.numberOfRows());
    ASSERT_EQ(4u, subMatrix.numberOfColumns());

    EXPECT_EQ(0, subMatrix(0,0));
    EXPECT_EQ(1, subMatrix(0,1));
    EXPECT_EQ(1, subMatrix(0,2));
    EXPECT_EQ(0, subMatrix(0,3));
    EXPECT_EQ(0, subMatrix(1,0));
    EXPECT_EQ(3, subMatrix(1,1));
    EXPECT_EQ(0, subMatrix(1,2));
    EXPECT_EQ(0, subMatrix(1,3));
    EXPECT_EQ(0, subMatrix(2,0));
    EXPECT_EQ(1, subMatrix(2,1));
    EXPECT_EQ(1, subMatrix(2,2));
    EXPECT_EQ(0, subMatrix(2,3));
}

template<class Matrix>
void MatricesGTest::testAssign() {
    // 1  0  1  0
    // 2  2  0  0
    // 3  0  0 -1
    // 0  0  1  4
    std::vector<Triplet> triplets = {{0,0,1}, {0,2,1}, {1,0,2}, {1,1,2}, {2,0,3}, {2,3,-1}, {3,2,1}, {3,3,4}};

    Matrix mat(4,4,triplets);

    // -1 1
    //  3 0
    std::vector<Triplet> subTriplets = {{0,0,-1}, {0,1,1}, {1,0,3}};
    Matrix sourceMat(2,2,subTriplets);

    std::vector<index> rows = {2,3};
    std::vector<index> columns = {0, 2};
    mat.assign(rows, columns, sourceMat);

    EXPECT_EQ(1, mat(0,0));
    EXPECT_EQ(0, mat(0,1));
    EXPECT_EQ(1, mat(0,2));
    EXPECT_EQ(0, mat(0,3));
    EXPECT_EQ(2, mat(1,0));
    EXPECT_EQ(2, mat(1,1));
    EXPECT_EQ(0, mat(1,2));
    EXPECT_EQ(0, mat(1,3));
    EXPECT_EQ(-1, mat(2,0));
    EXPECT_EQ(0, mat(2,1));
    EXPECT_EQ(1, mat(2,2));
    EXPECT_EQ(-1, mat(2,3));
    EXPECT_EQ(3, mat(3,0));
    EXPECT_EQ(0, mat(3,1));
    EXPECT_EQ(0, mat(3,2));
    EXPECT_EQ(4, mat(3,3));

    rows = {2,3};
    columns = {2,3};
    mat.assign(rows, columns, sourceMat);

    EXPECT_EQ(1, mat(0,0));
    EXPECT_EQ(0, mat(0,1));
    EXPECT_EQ(1, mat(0,2));
    EXPECT_EQ(0, mat(0,3));
    EXPECT_EQ(2, mat(1,0));
    EXPECT_EQ(2, mat(1,1));
    EXPECT_EQ(0, mat(1,2));
    EXPECT_EQ(0, mat(1,3));
    EXPECT_EQ(-1, mat(2,0));
    EXPECT_EQ(0, mat(2,1));
    EXPECT_EQ(-1, mat(2,2));
    EXPECT_EQ(1, mat(2,3));
    EXPECT_EQ(3, mat(3,0));
    EXPECT_EQ(0, mat(3,1));
    EXPECT_EQ(3, mat(3,2));
    EXPECT_EQ(0, mat(3,3));
}

template<class Matrix>
void MatricesGTest::testApply() {
    Matrix mat = Matrix::adjacencyMatrix(graph);

    mat.apply([&](double value) {return 2 * value;});
    graph.forEdges([&](index i, index j, double value) {
        EXPECT_EQ(2 * value, mat(i,j));
    });
}

template<class Matrix>
void MatricesGTest::testMatrixAddition() {
    std::vector<Triplet> triplets1;
    std::vector<Triplet> triplets2;

    for (index i = 0; i < 100; ++i) {
        triplets1.push_back({i,i,1});
        triplets2.push_back({i,i,(double)i});
    }

    triplets1.push_back({2,71, 1.8});
    triplets2.push_back({42,43,3.14});

    Matrix mat1(100, triplets1);
    Matrix mat2(100, triplets2);

    Matrix result = mat1 + mat2;
    ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

    EXPECT_EQ(0.0, result(10, 13));


    for (index i = 0; i < result.numberOfRows(); ++i) {
        EXPECT_EQ((i + 1), result(i, i));
    }
    EXPECT_EQ(1.8, result(2, 71));
    EXPECT_EQ(3.14, result(42, 43));

    EXPECT_EQ(0.0, result(3, 14));


    // rectangular matrix
    // n x m (n < m)
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    std::vector<Triplet> triplets = {{0,0,1.0}, {1,2,3.0}};


    mat1 = Matrix(2,5,triplets);

    // 0 0 1 0 0
    // 0 0 1 0 0
    triplets.clear();
    triplets = {{0,2,1.0}, {1,2,1.0}};

    mat2 = Matrix(2,5,triplets);

    result = mat1 + mat2;

    ASSERT_EQ(2u, result.numberOfRows());
    ASSERT_EQ(5u, result.numberOfColumns());

    EXPECT_EQ(1, result(0,0));
    EXPECT_EQ(1, result(0,2));
    EXPECT_EQ(4, result(1,2));

    EXPECT_EQ(0, result(0,1));
    EXPECT_EQ(0, result(1,4));

    // rectangular matrix
    // n x m (n > m)
    //
    // 1 0
    // 0 0
    // 0 3
    // 0 0
    // 0 0

    triplets.clear();
    triplets = {{0,0,1.0}, {2,1,3.0}};

    mat1 = Matrix(5,2,triplets);

    // 0 0
    // 0 0
    // 1 1
    // 0 0
    // 0 0
    triplets.clear();
    triplets = {{2,0,1.0}, {2,1,1.0}};

    mat2 = Matrix(5,2,triplets);

    result = mat1 + mat2;

    ASSERT_EQ(5u, result.numberOfRows());
    ASSERT_EQ(2u, result.numberOfColumns());

    EXPECT_EQ(1, result(0,0));
    EXPECT_EQ(1, result(2,0));
    EXPECT_EQ(4, result(2,1));

    EXPECT_EQ(0, result(0,1));
    EXPECT_EQ(0, result(4,1));
}

template<class Matrix>
void MatricesGTest::testMatrixSubtraction() {
    std::vector<Triplet> triplets1;
    std::vector<Triplet> triplets2;

    for (index i = 0; i < 100; ++i) {
        triplets1.push_back({i,i,1});
        triplets2.push_back({i,i,(double)i});
    }

    triplets1.push_back({2,71, 1.8});
    triplets2.push_back({42,43,3.14});

    Matrix mat1(100, triplets1);
    Matrix mat2(100, triplets2);

    Matrix result = mat2 - mat1;
    ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

    EXPECT_EQ(0.0, result(10, 13));


    for (index i = 0; i < result.numberOfRows(); ++i) {
        EXPECT_EQ(((int) i - 1), result(i, i));
    }
    EXPECT_EQ(-1.8, result(2, 71));
    EXPECT_EQ(3.14, result(42, 43));

    EXPECT_EQ(0.0, result(3, 14));


    // rectangular matrix
    // n x m (n < m)
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    std::vector<Triplet> triplets = {{0,0,1.0}, {1,2,3.0}};


    mat1 = Matrix(2,5,triplets);

    // 0 0 1 0 0
    // 0 0 1 0 0
    triplets.clear();
    triplets = {{0,2,1.0}, {1,2,1.0}};

    mat2 = Matrix(2,5,triplets);

    result = mat1 - mat2;

    ASSERT_EQ(2u, result.numberOfRows());
    ASSERT_EQ(5u, result.numberOfColumns());

    EXPECT_EQ(1, result(0,0));
    EXPECT_EQ(-1, result(0,2));
    EXPECT_EQ(2, result(1,2));

    EXPECT_EQ(0, result(0,1));
    EXPECT_EQ(0, result(1,4));

    // rectangular matrix
    // n x m (n > m)
    //
    // 1 0
    // 0 0
    // 0 3
    // 0 0
    // 0 0

    triplets.clear();
    triplets = {{0,0,1.0}, {2,1,3.0}};

    mat1 = Matrix(5,2,triplets);

    // 0 0
    // 0 0
    // 1 1
    // 0 0
    // 0 0
    triplets.clear();
    triplets = {{2,0,1.0}, {2,1,1.0}};

    mat2 = Matrix(5,2,triplets);

    result = mat1 - mat2;

    ASSERT_EQ(5u, result.numberOfRows());
    ASSERT_EQ(2u, result.numberOfColumns());

    EXPECT_EQ(1, result(0,0));
    EXPECT_EQ(-1, result(2,0));
    EXPECT_EQ(2, result(2,1));

    EXPECT_EQ(0, result(0,1));
    EXPECT_EQ(0, result(4,1));
}

template<class Matrix>
void MatricesGTest::testScalarMultiplication() {
    std::vector<Triplet> triplets;

    for (index i = 0; i < 100; ++i) {
        triplets.push_back({i,i,(double) i});
    }

    triplets.push_back({42,43,42.0});

    Matrix mat(100, triplets);
    mat *= 2;
    ASSERT_EQ(100u, mat.numberOfRows());
    ASSERT_EQ(100u, mat.numberOfColumns());

    for (index i = 0; i < 100; ++i) {
        EXPECT_EQ(i*2, mat(i, i));
    }
    EXPECT_EQ(84.0, mat(42, 43));
    EXPECT_EQ(0.0, mat(55, 99));

    mat *= 0.5;

    for (index i = 0; i < 100; ++i) {
        EXPECT_EQ(i, mat(i, i));
    }
    EXPECT_EQ(42.0, mat(42, 43));
    EXPECT_EQ(0.0, mat(55, 99));

    // rectangular matrix
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    triplets = {{0,0,1.0}, {1,2,3.0}};
    mat = Matrix(2,5,triplets);

    mat *= 2;

    EXPECT_EQ(2, mat(0,0));
    EXPECT_EQ(6, mat(1,2));
}

template<class Matrix>
void MatricesGTest::testMatrixDivisionOperator() {
    std::vector<Triplet> triplets;

    for (index i = 0; i < 100; ++i) {
        triplets.push_back({i,i, (double) i});
    }

    triplets.push_back({42,43,42.0});

    Matrix mat(100, triplets);
    mat /= (1.0 / 2.0);
    ASSERT_EQ(100u, mat.numberOfRows());
    ASSERT_EQ(100u, mat.numberOfColumns());

    for (index i = 0; i < 100; ++i) {
        EXPECT_EQ(i*2, mat(i, i));
    }
    EXPECT_EQ(84.0, mat(42, 43));
    EXPECT_EQ(0.0, mat(55, 99));

    mat /= 2;

    for (index i = 0; i < 100; ++i) {
        EXPECT_EQ(i, mat(i, i));
    }
    EXPECT_EQ(42.0, mat(42, 43));
    EXPECT_EQ(0.0, mat(55, 99));

    // rectangular matrix
    //
    // 1 0 0 0 0
    // 0 0 3 0 0
    triplets = {{0,0,1.0}, {1,2,3.0}};
    mat = Matrix(2,5,triplets);

    mat /= 2;

    EXPECT_EQ(0.5, mat(0,0));
    EXPECT_EQ(1.5, mat(1,2));
}

template<class Matrix>
void MatricesGTest::testMatrixVectorProduct() {
    std::vector<Triplet> triplets;

    for (index i = 0; i < 1000; ++i) {
        triplets.push_back({i,i, (double) i});
    }

    triplets.push_back({42,43,42.0});

    Vector vector(1000, 1.0);
    vector[500] = 3.5;

    Matrix mat(1000, triplets);

    Vector result = mat * vector;
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

    Matrix mat2(4, triplets);

    Vector v({1,2,3,0});
    Vector res = mat2 * v;
    ASSERT_EQ(mat2.numberOfRows(), res.getDimension());

    EXPECT_EQ(14, res[0]);
    EXPECT_EQ(6, res[1]);
    EXPECT_EQ(12, res[2]);
    EXPECT_EQ(-3, res[3]);

    // rectangular matrix
    //
    // 1 0 0 0 0
    // 0 0 3 0 0

    triplets = {{0,0,1}, {1,2,3}};
    mat = Matrix(2,5,triplets);

    v = {0,1,2,3,0};
    res = mat * v;

    ASSERT_EQ(2u, res.getDimension());
    EXPECT_EQ(0, res[0]);
    EXPECT_EQ(6, res[1]);
}

template<class Matrix>
void MatricesGTest::testMatrixMultiplication() {
    std::vector<Triplet> triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,2,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};

    //
    //				 1  2  3  0
    // 				 2  2  0  0
    // mat1 = mat2 = 3  0  3 -1
    //				 0  0 -1  4
    //
    Matrix mat1(4, triplets);
    ASSERT_EQ(4u, mat1.numberOfRows());
    ASSERT_EQ(4u, mat1.numberOfColumns());

    Matrix mat2(4, triplets);
    ASSERT_EQ(4u, mat2.numberOfRows());
    ASSERT_EQ(4u, mat2.numberOfColumns());

    //
    //			14  6  12  -3
    //			 6  8   6   0
    // result = 12  6  19  -7
    //			-3  0  -7  17
    //
    Matrix result = mat1 * mat2;
    ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

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


    // rectangular matrices
    //
    // 1 0 0 2
    // 0 0 1 0
    // 0 2 0 4
    triplets = {{0,0,1}, {0,3,2}, {1,2,1}, {2,1,2}, {2,3,4}};
    mat1 = Matrix(3,4,triplets);

    //
    // 1  0
    // 0  0
    // 0  0.5
    // 42 1

    triplets = {{0,0,1}, {2,1,0.5}, {3,0,42}, {3,1,1}};
    mat2 = Matrix(4,2, triplets);

    result = mat1 * mat2;
    ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat2.numberOfColumns(), result.numberOfColumns());

    EXPECT_EQ(85, result(0,0));
    EXPECT_EQ(2, result(0,1));
    EXPECT_EQ(0, result(1,0));
    EXPECT_EQ(0.5, result(1,1));
    EXPECT_EQ(168, result(2,0));
    EXPECT_EQ(4, result(2,1));
}

template<class Matrix>
void MatricesGTest::testBigMatrixMultiplication() {
    Matrix mat = Matrix::adjacencyMatrix(graph);
    Matrix result = mat * mat;
    ASSERT_EQ(mat.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat.numberOfColumns(), result.numberOfColumns());
}

template<class Matrix>
void MatricesGTest::testAdjacencyMatrix() {
    Graph G(6);
    G.addEdge(0,0);
    G.addEdge(0,1);
    G.addEdge(0,4);
    G.addEdge(1,2);
    G.addEdge(1,4);
    G.addEdge(2,3);
    G.addEdge(3,4);
    G.addEdge(3,5);

    Matrix mat = Matrix::adjacencyMatrix(G);

    // first row
    EXPECT_EQ(1, mat(0,0));
    EXPECT_EQ(1, mat(0,1));
    EXPECT_EQ(0, mat(0,2));
    EXPECT_EQ(0, mat(0,3));
    EXPECT_EQ(1, mat(0,4));
    EXPECT_EQ(0, mat(0,5));

    // third row
    EXPECT_EQ(0, mat(2,0));
    EXPECT_EQ(1, mat(2,1));
    EXPECT_EQ(0, mat(2,2));
    EXPECT_EQ(1, mat(2,3));
    EXPECT_EQ(0, mat(2,4));
    EXPECT_EQ(0, mat(2,5));

    // fifth row
    EXPECT_EQ(1, mat(4,0));
    EXPECT_EQ(1, mat(4,1));
    EXPECT_EQ(0, mat(4,2));
    EXPECT_EQ(1, mat(4,3));
    EXPECT_EQ(0, mat(4,4));
    EXPECT_EQ(0, mat(4,5));


    // directed, weighted G
    Graph dGraph(4, true, true);
    dGraph.addEdge(0,1,2);
    dGraph.addEdge(0,0, 42);
    dGraph.addEdge(2,3,-3);
    dGraph.addEdge(3,2,5);

    mat = Matrix::adjacencyMatrix(dGraph);
    ASSERT_EQ(dGraph.numberOfNodes(), mat.numberOfRows());
    ASSERT_EQ(dGraph.numberOfNodes(), mat.numberOfColumns());

    EXPECT_EQ(2, mat(0,1));
    EXPECT_EQ(0, mat(1,0));
    EXPECT_EQ(42, mat(0,0));
    EXPECT_EQ(-3, mat(2,3));
    EXPECT_EQ(5, mat(3,2));

    // read lesmis G
    METISGraphReader graphReader;
    G = graphReader.read("input/lesmis.graph");

    // create AdjacencyMatrix
    mat = Matrix::adjacencyMatrix(G);

    G.forNodes([&](node u) {
        G.forNodes([&](node v) {
            if (G.hasEdge(u,v)) {
                EXPECT_EQ(G.weight(u,v), mat(u,v));
            } else {
                EXPECT_EQ(0.0, mat(u,v));
            }
        });
    });
}

template<class Matrix>
void MatricesGTest::testDiagonalMatrix() {
    Vector diagonal = {1,0,4,-1};
    Matrix mat = Matrix::diagonalMatrix(diagonal);
    EXPECT_EQ(4u, mat.numberOfRows());
    EXPECT_EQ(4u, mat.numberOfColumns());

    EXPECT_EQ(1, mat(0,0));
    EXPECT_EQ(0, mat(1,1));
    EXPECT_EQ(4, mat(2,2));
    EXPECT_EQ(-1, mat(3,3));

    for (index i = 0; i < mat.numberOfRows(); ++i) {
        for (index j = 0; j < mat.numberOfColumns(); ++j) {
            if (i != j) {
                EXPECT_EQ(0, mat(i,j));
            }
        }
    }
}

template<class Matrix>
void MatricesGTest::testIncidenceMatrix() {
    Graph G = Graph(5, true);
    G.addEdge(0,1, 4.0);
    G.addEdge(0,2, 9.0);
    G.addEdge(0,3, 16.0);
    G.addEdge(2,3, 1.0);
    G.addEdge(4,1, 25.0);
    G.addEdge(4,4, 1.0);

    G.indexEdges();

    Matrix mat = Matrix::incidenceMatrix(G);
    ASSERT_EQ(G.numberOfNodes(), mat.numberOfRows());
    ASSERT_EQ(G.numberOfEdges(), mat.numberOfColumns());

    EXPECT_EQ(sqrt(G.weight(0,1)), mat(0,0));
    EXPECT_EQ(-sqrt(G.weight(0,1)), mat(1,0));
    for (uint64_t i = 2; i < mat.numberOfRows(); ++i) {
        EXPECT_EQ(0.0, mat(i, 0));
    }

    EXPECT_EQ(-sqrt(G.weight(0,2)), mat(2,1));

    EXPECT_EQ(-sqrt(G.weight(0,3)), mat(3,2));
    EXPECT_EQ(-sqrt(G.weight(2,3)), mat(3,3));

    for (uint64_t i = 0; i < mat.numberOfRows(); ++i) {
        EXPECT_EQ(0.0, mat(i, 5));
    }

    Vector row0 = mat.row(0);
    ASSERT_EQ(row0.getDimension(), mat.numberOfColumns());

    EXPECT_EQ(sqrt(G.weight(0,1)), row0[0]);
    EXPECT_EQ(sqrt(G.weight(0,2)), row0[1]);
    EXPECT_EQ(sqrt(G.weight(0,3)), row0[2]);
    for (uint64_t j = 3; j < row0.getDimension(); ++j) {
        EXPECT_EQ(0.0, row0[j]);
    }

    for (uint64_t j = 0; j < 5; ++j) {
        Vector column = mat.column(j);
        ASSERT_EQ(column.getDimension(), mat.numberOfRows());

        double sum = 0.0;
        for (uint64_t i = 0; i < column.getDimension(); ++i) {
            sum += column[i];
        }

        EXPECT_EQ(0.0, sum);
    }

    Vector column5 = mat.column(5);
    ASSERT_EQ(column5.getDimension(), mat.numberOfRows());

    for (uint64_t i = 0; i < column5.getDimension(); ++i) {
        EXPECT_EQ(0.0, column5[i]);
    }

    Vector v = {12, 3, 9, 28, 0, -1};

    Vector result = mat * v;
    ASSERT_EQ(result.getDimension(), mat.numberOfRows());

    EXPECT_EQ(69, result[0]);
    EXPECT_EQ(-24, result[1]);
    EXPECT_EQ(19, result[2]);
    EXPECT_EQ(-64, result[3]);
    EXPECT_EQ(0, result[4]);
}


template<class Matrix>
void MatricesGTest::testLaplacianOfGraph() {
    METISGraphReader graphReader;
    Matrix mat = Matrix::laplacianMatrix(graph);
    EXPECT_TRUE(MatrixTools::isLaplacian(mat));

    mat = Matrix::laplacianMatrix(graphReader.read("input/power.graph"));
    EXPECT_TRUE(MatrixTools::isLaplacian(mat));
}



TEST_F(MatricesGTest, testDimension) {
    testDimension<DynamicMatrix>();
    testDimension<CSRMatrix>();
    testDimension<DenseMatrix>();
}

TEST_F(MatricesGTest, testNNZInRow) {
    testNNZInRow<DynamicMatrix>();
    testNNZInRow<CSRMatrix>();
    testNNZInRow<DenseMatrix>();
}

TEST_F(MatricesGTest, testRowAndColumnAccess) {
    testRowAndColumnAccess<DynamicMatrix>();
    testRowAndColumnAccess<CSRMatrix>();
    testRowAndColumnAccess<DenseMatrix>();
}

TEST_F(MatricesGTest, testDiagonalVector) {
    testDiagonalVector<DynamicMatrix>();
    testDiagonalVector<CSRMatrix>();
    testDiagonalVector<DenseMatrix>();
}

TEST_F(MatricesGTest, testTranspose) {
    testTranspose<DynamicMatrix>();
    testTranspose<CSRMatrix>();
    testTranspose<DenseMatrix>();
}

TEST_F(MatricesGTest, testExtract) {
    testExtract<DynamicMatrix>();
    testExtract<CSRMatrix>();
}

TEST_F(MatricesGTest, testAssign) {
    testAssign<DynamicMatrix>();
    testAssign<CSRMatrix>();
    testAssign<DenseMatrix>();
}

TEST_F(MatricesGTest, testApply) {
    testApply<DynamicMatrix>();
    testApply<CSRMatrix>();
}

TEST_F(MatricesGTest, testMatrixAddition) {
    testMatrixAddition<DynamicMatrix>();
    testMatrixAddition<CSRMatrix>();
    testMatrixAddition<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixSubtraction) {
    testMatrixSubtraction<DynamicMatrix>();
    testMatrixSubtraction<CSRMatrix>();
    testMatrixSubtraction<DenseMatrix>();
}

TEST_F(MatricesGTest, testScalarMultiplication) {
    testScalarMultiplication<DynamicMatrix>();
    testScalarMultiplication<CSRMatrix>();
    testScalarMultiplication<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixDivisionOperator) {
    testMatrixDivisionOperator<DynamicMatrix>();
    testMatrixDivisionOperator<CSRMatrix>();
    testMatrixDivisionOperator<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixVectorProduct) {
    testMatrixVectorProduct<DynamicMatrix>();
    testMatrixVectorProduct<CSRMatrix>();
    testMatrixVectorProduct<DenseMatrix>();
}

TEST_F(MatricesGTest, testMatrixMultiplication) {
    testMatrixMultiplication<DynamicMatrix>();
    testMatrixMultiplication<CSRMatrix>();
    testMatrixMultiplication<DenseMatrix>();
}

TEST_F(MatricesGTest, testBigMatrixMultiplcation) {
    testBigMatrixMultiplication<DynamicMatrix>();
    testBigMatrixMultiplication<CSRMatrix>();
}

TEST_F(MatricesGTest, testAdjacencyMatrixOfGraph) {
    testAdjacencyMatrix<DynamicMatrix>();
    testAdjacencyMatrix<CSRMatrix>();
}

TEST_F(MatricesGTest, testDiagonalMatrix) {
    testDiagonalMatrix<DynamicMatrix>();
    testDiagonalMatrix<CSRMatrix>();
}

TEST_F(MatricesGTest, testIncidenceMatrix) {
    testIncidenceMatrix<DynamicMatrix>();
    testIncidenceMatrix<CSRMatrix>();
}

TEST_F(MatricesGTest, testLaplacianMatrixOfGraph) {
    testLaplacianOfGraph<DynamicMatrix>();
    testLaplacianOfGraph<CSRMatrix>();
}

} /* namespace NetworKit */
