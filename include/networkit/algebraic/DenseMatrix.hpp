// no-networkit-format
/*
 * DenseMatrix.hpp
 *
 *  Created on: Nov 25, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_ALGEBRAIC_DENSE_MATRIX_HPP_
#define NETWORKIT_ALGEBRAIC_DENSE_MATRIX_HPP_

#include <cassert>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>
#include <networkit/algebraic/Vector.hpp>

namespace NetworKit {

/**
 * @ingroup algebraic
 * Represents a dense matrix. Use this matrix to run LU decompositions and LU solves.
 * Note that most matrices are rather sparse s.t. CSRMatrix might be a better representation.
 */
class DenseMatrix final {
private:
    count nRows;
    count nCols;
    std::vector<double> entries;
    double zero;

public:
    /** Default constructor */
    DenseMatrix();

    /**
     * Constructs the DenseMatrix with size @a dimension x @a dimension.
     * @param dimension Defines how many rows and columns this matrix has.
     * @param zero The zero element (default is 0.0).
     */
    DenseMatrix(count dimension, double zero = 0.0);

    /**
     * Constructs the DenseMatrix with size @a nRows x @a nCols.
     * @param nRows Number of rows.
     * @param nCols Number of columns.
     * @param zero The zero element (default is 0.0).
     */
    DenseMatrix(count nRows, count nCols, double zero = 0.0);

    /**
     * Constructs the @a dimension x @a dimension DenseMatrix from the elements at position @a positions with values @values.
     * @param dimension Defines how many rows and columns this matrix has.
     * @param triplets The nonzero elements.
     * @param zero The zero element (default is 0.0).
     */
    DenseMatrix(count dimension, const std::vector<Triplet>& triplets, double zero = 0.0);

    /**
     * Constructs the @a nRows x @a nCols DenseMatrix from the elements at position @a positions with values @values.
     * @param nRows Defines how many rows this matrix has.
     * @param nCols Defines how many columns this matrix has.
     * @param triplets The nonzero elements.
     * @param zero The zero element (default is 0.0).
     */
    DenseMatrix(count nRows,count nCols, const std::vector<Triplet>& triplets, double zero = 0.0);

    /**
     * Constructs an instance of DenseMatrix given the number of rows (@a nRows) and the number of columns (@a nCols) and its
     * values (@a entries).
     * @param nRows Number of rows.
     * @param nCols Number of columns.
     * @param entries Entries of the matrix.
     * @param zero The zero element (default is 0.0).
     * @note The size of the @a entries vector should be equal to @a nRows * @a nCols.
     */
    DenseMatrix(count nRows, count nCols, const std::vector<double>& entries, double zero = 0.0);

    /** Default destructor */
    ~DenseMatrix() = default;

    /** Default copy constructor */
    DenseMatrix (const DenseMatrix &other) = default;

    /** Default move constructor */
    DenseMatrix (DenseMatrix &&other) = default;

    /** Default copy assignment operator */
    DenseMatrix& operator=(DenseMatrix &&other) = default;

    /** Default move assignment operator */
    DenseMatrix& operator=(const DenseMatrix &other) = default;

    /**
     * @return Number of rows.
     */
    inline count numberOfRows() const {
        return nRows;
    }

    /**
     * @return Number of columns.
     */
    inline count numberOfColumns() const {
        return nCols;
    }

    /**
     * Returns the zero element of the matrix.
     */
    inline double getZero() const {
        return zero;
    }

    /**
     * @param i The row index.
     * @return Number of non-zeros in row @a i.
     * @note This function is linear in the number of columns of the matrix.
     */
    count nnzInRow(index i) const;

    /**
     * @return Number of non-zeros in this matrix.
     * @note This function takes nRows * nCols operations.
     */
    count nnz() const;

    /**
     * @return Value at matrix position (i,j).
     */
    double operator()(index i, index j) const;

    /**
     * Set the matrix at position (@a i, @a j) to @a value.
     */
    void setValue(index i, index j, double value);


    /**
     * @return Row @a i of this matrix as vector.
     */
    Vector row(index i) const;

    /**
     * @return Column @a j of this matrix as vector.
     */
    Vector column(index j) const;

    /**
     * @return The main diagonal of this matrix.
     */
    Vector diagonal() const;

    /**
     * Adds this matrix to @a other and returns the result.
     * @return The sum of this matrix and @a other.
     */
    DenseMatrix operator+(const DenseMatrix &other) const;

    /**
     * Adds @a other to this matrix.
     * @return Reference to this matrix.
     */
    DenseMatrix& operator+=(const DenseMatrix &other);

    /**
     * Subtracts @a other from this matrix and returns the result.
     * @return The difference of this matrix and @a other.
     *
     */
    DenseMatrix operator-(const DenseMatrix &other) const;

    /**
     * Subtracts @a other from this matrix.
     * @return Reference to this matrix.
     */
    DenseMatrix& operator-=(const DenseMatrix &other);

    /**
     * Multiplies this matrix with a scalar specified in @a scalar and returns the result.
     * @return The result of multiplying this matrix with @a scalar.
     */
    DenseMatrix operator*(double scalar) const;

    /**
     * Multiplies this matrix with a scalar specified in @a scalar.
     * @return Reference to this matrix.
     */
    DenseMatrix& operator*=(double scalar);

    /**
     * Multiplies this matrix with @a vector and returns the result.
     * @return The result of multiplying this matrix with @a vector.
     */
    Vector operator*(const Vector &vector) const;

    /**
     * Multiplies this matrix with @a other and returns the result in a new matrix.
     * @return The result of multiplying this matrix with @a other.
     */
    DenseMatrix operator*(const DenseMatrix &other) const;

    /**
     * Divides this matrix by a divisor specified in @a divisor and returns the result in a new matrix.
     * @return The result of dividing this matrix by @a divisor.
     */
    DenseMatrix operator/(double divisor) const;

    /**
     * Divides this matrix by a divisor specified in @a divisor.
     * @return Reference to this matrix.
     */
    DenseMatrix& operator/=(double divisor);

    /**
     * Transposes this matrix and returns it.
     */
    DenseMatrix transpose() const;

    /**
     * Extracts a matrix with rows and columns specified by @a rowIndices and @a columnIndices from this matrix.
     * The order of rows and columns is equal to the order in @a rowIndices and @a columnIndices. It is also
     * possible to specify a row or column more than once to get duplicates.
     * @param rowIndices
     * @param columnIndices
     */
    DenseMatrix extract(const std::vector<index>& rowIndices, const std::vector<index>& columnIndices) const;

    /**
     * Assign the contents of the matrix @a source to this matrix at rows and columns specified by @a rowIndices and
     * @a columnIndices. That is, entry (i,j) of @a source is assigned to entry (rowIndices[i], columnIndices[j]) of
     * this matrix. Note that the dimensions of @rowIndices and @a columnIndices must coincide with the number of rows
     * and columns of @a source.
     * @param rowIndices
     * @param columnIndices
     * @param source
     */
    void assign(const std::vector<index>& rowIndices, const std::vector<index>& columnIndices, const DenseMatrix& source);

    /**
     * Applies the unary function @a unaryElementFunction to each value in the matrix. Note that it must hold that the
     * function applied to the zero element of this matrix returns the zero element.
     * @param unaryElementFunction
     */
    template<typename F>
    void apply(F unaryElementFunction);

    /**
     * Decomposes the given @a matrix into lower L and upper U matrix (in-place).
     * @param matrix The matrix to decompose into LU.
     */
    static void LUDecomposition(DenseMatrix &matrix);

    /**
     * Computes the solution vector x to the system @a LU * x = @a b where @a LU is a matrix decomposed into L and U.
     * @param LU Matrix decomposed into lower L and upper U matrix.
     * @param b Right-hand side.
     * @return Solution vector x to the linear equation system LU * x = b.
     */
    static Vector LUSolve(const DenseMatrix &LU, const Vector &b);

    /**
     * Computes @a A @a binaryOp @a B on the elements of matrix @a A and matrix @a B.
     * @param A
     * @param B
     * @param binaryOp Function handling (double, double) -> double
     * @return @a A @a binaryOp @a B.
     * @note @a A and @a B must have the same dimensions.
     */
    template<typename L> static DenseMatrix binaryOperator(const DenseMatrix &A, const DenseMatrix &B, L binaryOp);

    /**
     * Iterate over all non-zero elements of row @a row in the matrix and call handler(index column, double value)
     */
    template<typename L> void forElementsInRow(index row, L handle) const;

    /**
     * Iterate in parallel over all non-zero elements of row @a row in the matrix and call handler(index column, double value)
     */
    template<typename L> void parallelForElementsInRow(index row, L handle) const;

    /**
     * Iterate over all non-zero elements of the matrix in row order and call handler (lambda closure).
     */
    template<typename L> void forElementsInRowOrder(L handle) const;

    /**
     * Iterate in parallel over all rows and call handler (lambda closure) on non-zero elements of the matrix.
     */
    template<typename L> void parallelForElementsInRowOrder(L handle) const;

    /**
     * Iterate over all non-zero elements of row @a row in the matrix and call handler(index column, double value).
     * @note This is a DenseMatrix! Therefore this operation needs O(numberOfRows()) time regardless of the number of
     * non-zeros actually present.
     */
    template<typename L> void forNonZeroElementsInRow(index row, L handle) const;

    /**
     * Iterate in parallel over all non-zero elements of row @a row in the matrix and call handler(index column, double value)
     * @note This is a DenseMatrix! Therefore this operation needs O(numberOfRows()) sequential time regardless of the number
     * of non-zeros actually present.
     */
    template<typename L> void parallelForNonZeroElementsInRow(index row, L handle) const;

    /**
     * Iterate over all non-zero elements of the matrix in row order and call handler (lambda closure).
     * @note This is a DenseMatrix! Therefore this operation needs O(numberOfRows() * numberOfColumns()) time regardless of the
     *  number of non-zeros actually present.
     */
    template<typename L> void forNonZeroElementsInRowOrder(L handle) const;

    /**
     * Iterate in parallel over all rows and call handler (lambda closure) on non-zero elements of the matrix.
     * @note This is a DenseMatrix! Therefore this operation needs O(numberOfRows() * numberOfColumns()) sequential time regardless
     *  of the number of non-zeros actually present.
     */
    template<typename L> void parallelForNonZeroElementsInRowOrder(L handle) const;
};

template<typename F>
void DenseMatrix::apply(F unaryElementFunction) {
#pragma omp parallel for
    for (omp_index k = 0; k < static_cast<omp_index>(entries.size()); ++k) {
        entries[k] = unaryElementFunction(entries[k]);
    }
}

template<typename L> inline DenseMatrix DenseMatrix::binaryOperator(const DenseMatrix &A, const DenseMatrix &B, L binaryOp) {
    assert(A.nRows == B.nRows && A.nCols == B.nCols);

    std::vector<double> resultEntries(A.numberOfRows() * A.numberOfColumns(), 0.0);

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(A.numberOfRows()); ++i) {
        index offset = i * A.numberOfColumns();
        for (index j = offset; j < offset + A.numberOfColumns(); ++j) {
            resultEntries[j] = binaryOp(A.entries[j], B.entries[j]);
        }
    }

    return DenseMatrix(A.numberOfRows(), A.numberOfColumns(), resultEntries);
}

template<typename L>
inline void DenseMatrix::forElementsInRow(index i, L handle) const {
    index offset = i * numberOfColumns();
    for (index k = offset, j = 0; k < offset + numberOfColumns(); ++k, ++j) {
        handle(j, entries[k]);
    }
}

template<typename L>
inline void DenseMatrix::parallelForElementsInRow(index i, L handle) const {
    index offset = i * numberOfColumns();
#pragma omp parallel for
    for (omp_index j = 0; j < static_cast<omp_index>(numberOfColumns()); ++j) {
        handle(j, entries[offset + j]);
    }
}

template<typename L>
inline void DenseMatrix::forElementsInRowOrder(L handle) const {
    for (index i = 0; i < nRows; ++i) {
        index offset = i * numberOfColumns();
        for (index k = offset, j = 0; k < offset + numberOfColumns(); ++k, ++j) {
            handle(i, j, entries[k]);
        }
    }
}

template<typename L>
inline void DenseMatrix::parallelForElementsInRowOrder(L handle) const {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
        index offset = i * numberOfColumns();
        for (index k = offset, j = 0; k < offset + numberOfColumns(); ++k, ++j) {
            handle(i, j, entries[k]);
        }
    }
}

template<typename L>
inline void DenseMatrix::forNonZeroElementsInRow(index row, L handle) const {
    for (index j = 0, k = row * numberOfColumns(); j < numberOfColumns(); ++j, ++k) {
        if (entries[k] != getZero()) {
            handle(j, entries[k]);
        }
    }
}

template<typename L>
inline void DenseMatrix::parallelForNonZeroElementsInRow(index row, L handle) const {
#pragma omp parallel for
    for (omp_index j = 0; j < static_cast<omp_index>(numberOfColumns()); ++j) {
        index k = row * numberOfColumns() + j;
        if (entries[k] != getZero()) {
            handle(j, entries[k]);
        }
    }
}

template<typename L>
inline void DenseMatrix::forNonZeroElementsInRowOrder(L handle) const {
    for (index i = 0; i < numberOfRows(); ++i) {
        for (index j = 0, k = i * numberOfColumns(); j < numberOfColumns(); ++j, ++k) {
            if (entries[k] != getZero()) {
                handle(i,j,entries[k]);
            }
        }
    }
}

template<typename L>
inline void DenseMatrix::parallelForNonZeroElementsInRowOrder(L handle) const {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(numberOfRows()); ++i) {
        for (index j = 0, k = i * numberOfColumns(); j < numberOfColumns(); ++j, ++k) {
            if (entries[k] != getZero()) {
                handle(i,j,entries[k]);
            }
        }
    }
}

} /* namespace NetworKit */

#endif // NETWORKIT_ALGEBRAIC_DENSE_MATRIX_HPP_
