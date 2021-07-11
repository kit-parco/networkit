// no-networkit-format
/*
 * GraphBLAS.h
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_ALGEBRAIC_GRAPH_BLAS_HPP_
#define NETWORKIT_ALGEBRAIC_GRAPH_BLAS_HPP_

#include <limits>
#include <networkit/algebraic/Semirings.hpp>
#include <networkit/algebraic/SparseAccumulator.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>
#include <networkit/algebraic/Vector.hpp>

/**
 * @ingroup algebraic
 * Implements the GraphBLAS interface. For more information visit https://graphblas.org.
 */
namespace GraphBLAS {

// ****************************************************
// 						Operations
// ****************************************************

/**
 * Computes binOp(A(i,j), B(i,j)) for all i,j element-wise. Note that the dimensions of
 * @a A and @a B must coincide and that the zero must be the same.
 * @param A
 * @param B
 * @param binOp
 * @return The resulting matrix.
 */
template<class SemiRing, class Matrix, typename L>
Matrix eWiseBinOp(const Matrix& A, const Matrix& B, L binOp) {
    assert(A.numberOfRows() == B.numberOfRows() && A.numberOfColumns() == B.numberOfColumns());
    assert(A.getZero() == B.getZero() && A.getZero() == SemiRing::zero());

    std::vector<int64_t> columnPointer(A.numberOfColumns(), -1);
    std::vector<double> Arow(A.numberOfColumns(), SemiRing::zero());
    std::vector<double> Brow(A.numberOfColumns(), SemiRing::zero());

    std::vector<NetworKit::Triplet> triplets;

    for (NetworKit::index i = 0; i < A.numberOfRows(); ++i) {
        NetworKit::index listHead = 0;
        NetworKit::count nnz = 0;

        // search for nonZeros in matrix A
        A.forNonZeroElementsInRow(i, [&](NetworKit::index j, double value) {
            Arow[j] = value;

            columnPointer[j] = listHead;
            listHead = j;
            nnz++;
        });

        // search for nonZeros in matrix B
        B.forNonZeroElementsInRow(i, [&](NetworKit::index j, double value) {
            Brow[j] = value;

            if (columnPointer[j] == -1) { // matrix A does not have a nonZero entry in column j
                columnPointer[j] = listHead;
                listHead = j;
                nnz++;
            }
        });

        // apply operator on the found nonZeros in A and B
        for (NetworKit::count k = 0; k < nnz; ++k) {
            double value = binOp(Arow[listHead], Brow[listHead]);
            if (value != SemiRing::zero()) {
                triplets.push_back({i,listHead,value});
            }

            NetworKit::index temp = listHead;
            listHead = columnPointer[listHead];

            // reset for next row
            columnPointer[temp] = -1;
            Arow[temp] = SemiRing::zero();
            Brow[temp] = SemiRing::zero();
        }

        nnz = 0;
    }

    return Matrix(A.numberOfRows(), A.numberOfColumns(), triplets, A.getZero());
}

/**
 * Computes the matrix-matrix multiplication of @a A and @a B. Note that
 * A.numberOfColumns() must be equal to B.numberOfRows() and the zero elements
 * must be the same. The default Semiring is the ArithmeticSemiring.
 * @param A
 * @param B
 * @return The result of the multiplication A * B.
 */
template<class SemiRing = ArithmeticSemiring, class Matrix>
Matrix MxM(const Matrix& A, const Matrix& B) {
    assert(A.numberOfColumns() == B.numberOfRows());
    assert(A.getZero() == SemiRing::zero() && B.getZero() == SemiRing::zero());

    std::vector<NetworKit::Triplet> triplets;
    NetworKit::SparseAccumulator spa(B.numberOfRows());
    for (NetworKit::index i = 0; i < A.numberOfRows(); ++i) {
        A.forNonZeroElementsInRow(i, [&](NetworKit::index k, double w1) {
            B.forNonZeroElementsInRow(k, [&](NetworKit::index j, double w2) {
                spa.scatter(SemiRing::mult(w1,w2), j, *SemiRing::add);
            });
        });

        spa.gather([&](NetworKit::index i, NetworKit::index j, double value){
            triplets.push_back({i,j,value});
        });

        spa.increaseRow();
    }

    return Matrix(A.numberOfRows(), B.numberOfColumns(), triplets, A.getZero());
}

/**
 * Computes the matrix-matrix multiplication of @a A and @a B and adds it to @a C where
 * the add operation is that of the specified Semiring (i.e. C(i,j) = SemiRing::add(C(i,j), (A*B)(i,j))).
 * The default Semiring is the ArithmeticSemiring.
 * @param A
 * @param B
 * @param C
 */
template<class SemiRing = ArithmeticSemiring, class Matrix>
void MxM(const Matrix& A, const Matrix& B, Matrix& C) {
    assert(A.numberOfColumns() == B.numberOfRows() && A.numberOfRows() == C.numberOfRows() && B.numberOfColumns() == C.numberOfColumns());
    assert(A.getZero() == SemiRing::zero() && B.getZero() == SemiRing::zero() && C.getZero() == SemiRing::zero());

    std::vector<NetworKit::Triplet> triplets;
    NetworKit::SparseAccumulator spa(B.numberOfRows());
    for (NetworKit::index i = 0; i < A.numberOfRows(); ++i) {
        A.forNonZeroElementsInRow(i, [&](NetworKit::index k, double w1) {
            B.forNonZeroElementsInRow(k, [&](NetworKit::index j, double w2) {
                spa.scatter(SemiRing::mult(w1,w2), j, *SemiRing::add);
            });
        });

        spa.gather([&](NetworKit::index i, NetworKit::index j, double value){
            triplets.push_back({i,j,value});
        });

        spa.increaseRow();
    }

    Matrix temp(A.numberOfRows(), B.numberOfRows(), triplets, A.getZero());
    C = eWiseBinOp<SemiRing, Matrix>(C, temp, *SemiRing::add);
}

/**
 * Computes the matrix-matrix multiplication of @a A and @a B and adds it to @a C where
 * the add operation is specified by the binary function @a accum (i.e. C(i,j) = accum(C(i,j), (A*B)(i,j))).
 * The default Semiring is the ArithmeticSemiring.
 * @param A
 * @param B
 * @param C
 * @param accum
 */
template<class SemiRing = ArithmeticSemiring, typename F, class Matrix>
void MxM(const Matrix& A, const Matrix& B, Matrix& C, F accum) {
    assert(A.numberOfColumns() == B.numberOfRows() && A.numberOfRows() == C.numberOfRows() && B.numberOfColumns() == C.numberOfColumns());
    assert(A.getZero() == SemiRing::zero() && B.getZero() == SemiRing::zero() && C.getZero() == SemiRing::zero());

    std::vector<NetworKit::Triplet> triplets;
    NetworKit::SparseAccumulator spa(B.numberOfRows());
    for (NetworKit::index i = 0; i < A.numberOfRows(); ++i) {
        A.forNonZeroElementsInRow(i, [&](NetworKit::index k, double w1) {
            B.forNonZeroElementsInRow(k, [&](NetworKit::index j, double w2) {
                spa.scatter(SemiRing::mult(w1,w2), j, *SemiRing::add);
            });
        });

        spa.gather([&](NetworKit::index i, NetworKit::index j, double value){
            triplets.push_back({i,j,value});
        });

        spa.increaseRow();
    }

    Matrix temp(A.numberOfRows(), B.numberOfRows(), triplets, A.getZero());
    C = eWiseBinOp<SemiRing, Matrix>(C, temp, accum);
}

/**
 * Computes the matrix-vector product of matrix @a A and Vector @a v. The default Semiring is the ArithmeticSemiring.
 * @param A
 * @param v
 */
template<class SemiRing = ArithmeticSemiring, class Matrix>
NetworKit::Vector MxV(const Matrix& A, const NetworKit::Vector& v) {
    assert(!v.isTransposed());
    assert(A.numberOfColumns() == v.getDimension());
    assert(A.getZero() == SemiRing::zero());
    NetworKit::Vector result(A.numberOfRows(), A.getZero());

    A.parallelForNonZeroElementsInRowOrder([&](NetworKit::index i, NetworKit::index j, double value) {
        result[i] = SemiRing::add(result[i], SemiRing::mult(value, v[j]));
    });

    return result;
}

/**
 * Computes the matrix-vector product of matrix @a A and Vector @a v and adds it to @a c where the add operation
 * is that of the specified Semiring (i.e. c[i] = SemiRing::add(c[i], (A*v)[i]). The default Semiring is the
 * ArithmeticSemiring.
 * @param A
 * @param v
 * @param c
 */
template<class SemiRing = ArithmeticSemiring, class Matrix>
void MxV(const Matrix& A, const NetworKit::Vector& v, NetworKit::Vector& c) {
    assert(!v.isTransposed());
    assert(A.numberOfColumns() == v.getDimension());
    assert(A.getZero() == SemiRing::zero());

    A.parallelForNonZeroElementsInRowOrder([&](NetworKit::index i, NetworKit::index j, double value) {
        c[i] = SemiRing::add(c[i], SemiRing::mult(value, v[j]));
    });
}

/**
 * Computes the matrix-vector product of matrix @a A and Vector @a v and adds it to @a c where the add operation
 * is that of the specified binary function @a accum (i.e. c[i] = accum(c[i], (A*v)[i]). The default Semiring is the
 * ArithmeticSemiring.
 * @param A
 * @param v
 * @param c
 */
template<class SemiRing = ArithmeticSemiring, typename F, class Matrix>
void MxV(const Matrix& A, const NetworKit::Vector& v, NetworKit::Vector& c, F accum) {
    assert(!v.isTransposed());
    assert(A.numberOfColumns() == v.getDimension());
    assert(A.getZero() == SemiRing::zero());

    A.parallelForNonZeroElementsInRowOrder([&](NetworKit::index i, NetworKit::index j, double value) {
        c[i] = accum(c[i], SemiRing::mult(value, v[j]));
    });
}

/**
 * Computes SemiRing::add(A(i,j), B(i,j)) for all i,j element-wise and returns the resulting matrix. The default
 * Semiring is the ArithmeticSemiring.
 * @param A
 * @param B
 */
template<class SemiRing = ArithmeticSemiring, class Matrix>
Matrix eWiseAdd(const Matrix& A, const Matrix& B) {
    return eWiseBinOp<SemiRing, Matrix>(A, B, [](const double a, const double b) {return SemiRing::add(a,b);});
}

/**
 * Computes SemiRing::mult(A(i,j), B(i,j)) for all i,j element-wise and returns the resulting matrix. The default
 * Semiring is the ArithmeticSemiring.
 * @param A
 * @param B
 * @return
 */
template<class SemiRing = ArithmeticSemiring, class Matrix>
Matrix eWiseMult(const Matrix& A, const Matrix& B) {
    return eWiseBinOp<SemiRing, Matrix>(A, B, [](const double a, const double b) {return SemiRing::mult(a,b);});
}

/**
 * Computes the row-reduction of the @a matrix and returns the result as a vector. That is, the elements of each row
 * are summed up to form the respective entry in the result vector. The add operator is that of the specified
 * Semiring. The default Semiring is the ArithmeticSemiring.
 * @param matrix
 */
template<class SemiRing = ArithmeticSemiring, class Matrix>
NetworKit::Vector rowReduce(const Matrix& matrix) {
    assert(matrix.getZero() == SemiRing::zero());
    NetworKit::Vector rowReduction(matrix.numberOfRows(), 0.0);

#pragma omp parallel for
    for (NetworKit::omp_index i = 0; i < static_cast<NetworKit::omp_index>(matrix.numberOfRows()); ++i) {
        matrix.forNonZeroElementsInRow(i, [&](NetworKit::index, double value) {
            rowReduction[i] = SemiRing::add(rowReduction[i], value);
        });
    }

    return rowReduction;
}

/**
 * Computes the column-reduction of the @a matrix and returns the result as a Vector. That is, the elements of each
 * column are summed up to form the respective entry in the result Vector. The add operator is that of the specified
 * Semiring. The default Semiring is the ArithmeticSemiring.
 * @param matrix
 */
template<class SemiRing = ArithmeticSemiring, class Matrix>
NetworKit::Vector columnReduce(const Matrix& matrix) {
    assert(matrix.getZero() == SemiRing::zero());
    NetworKit::Vector columnReduction(matrix.numberOfColumns(), 0.0);

    matrix.forNonZeroElementsInRowOrder([&](NetworKit::index, NetworKit::index j, double value) {
        columnReduction[j] = SemiRing::add(columnReduction[j], value);
    });

    return columnReduction;
}

}



#endif // NETWORKIT_ALGEBRAIC_GRAPH_BLAS_HPP_
