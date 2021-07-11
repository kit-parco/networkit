/*
 * CSRGeneralMatrix.hpp
 *
 *  Created on: May 6, 2015
 *     Authors: Michael Wegner <michael.wegner@student.kit.edu>
 *              Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_ALGEBRAIC_CSR_GENERAL_MATRIX_HPP_
#define NETWORKIT_ALGEBRAIC_CSR_GENERAL_MATRIX_HPP_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <omp.h>
#include <stdexcept>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/graph/Graph.hpp>

#include <tlx/unused.hpp>

namespace NetworKit {

/**
 * @ingroup algebraic
 * The CSRGeneralMatrix class represents a sparse matrix stored in CSR-Format
 * (i.e. compressed sparse row).
 * If speed is important, use this CSRGeneralMatrix instead of the Matrix class.
 */
template <class ValueType>
class CSRGeneralMatrix {
    std::vector<index> rowIdx, columnIdx;
    std::vector<ValueType> nonZeros;

    count nRows, nCols;
    bool isSorted;
    ValueType zero;

    /**
     * Quicksort algorithm on columnIdx between [@a left, @a right].
     * @param left
     * @param right
     */
    void quicksort(index left, index right) {
        if (left >= right)
            return;
        index pivotIdx = partition(left, right);
        if (pivotIdx != 0) {
            quicksort(left, pivotIdx - 1);
        }
        quicksort(pivotIdx + 1, right);
    }

    /**
     * Partitions columnIdx between [@a left, @a right] after selecting the pivot
     * in the middle.
     * @param left
     * @param right
     * @return The pivot.
     */
    index partition(index left, index right) {
        index mid = (left + right) / 2;
        index pivot = columnIdx[mid];
        std::swap(columnIdx[mid], columnIdx[right]);
        std::swap(nonZeros[mid], nonZeros[right]);

        index i = left;
        for (index j = left; j < right; ++j) {
            if (columnIdx[j] <= pivot) {
                std::swap(columnIdx[i], columnIdx[j]);
                std::swap(nonZeros[i], nonZeros[j]);
                ++i;
            }
        }
        std::swap(columnIdx[i], columnIdx[right]);
        std::swap(nonZeros[i], nonZeros[right]);
        return i;
    }

    /**
     * Binary search the sorted columnIdx vector between [@a left, @a right]
     * for column @a j.
     * If @a j is not present, the index that is immediately left of the place
     * where @a j would be is returned.
     * @param left
     * @param right
     * @param j
     * @return The position of column @a j in columnIdx or the element immediately
     * to the left of the place where @a j would be.
     */
    index binarySearchColumns(index left, index right, index j) const {
        assert(sorted());
        const auto it = std::lower_bound(columnIdx.begin() + left, columnIdx.begin() + right, j);
        if (it == columnIdx.end() || *it != j)
            return none;
        return it - columnIdx.begin();
    }

public:
    /** Default constructor */
    CSRGeneralMatrix()
        : rowIdx(0), columnIdx(0), nonZeros(0), nRows(0), nCols(0), isSorted(true), zero(0) {}

    /**
     * Constructs the CSRGeneralMatrix with size @a dimension x @a dimension.
     * @param dimension Defines how many rows and columns this matrix has.
     * @param zero The zero element (default = 0).
     */
    CSRGeneralMatrix(count dimension, ValueType zero = 0)
        : rowIdx(dimension + 1), columnIdx(0), nonZeros(0), nRows(dimension), nCols(dimension),
          isSorted(true), zero(zero) {}

    /**
     * Constructs the CSRGeneralMatrix with size @a nRows x @a nCols.
     * @param nRows Number of rows.
     * @param nCols Number of columns.
     * @param zero The zero element (default = 0).
     */
    CSRGeneralMatrix(count nRows, count nCols, ValueType zero = 0)
        : rowIdx(nRows + 1), columnIdx(0), nonZeros(0), nRows(nRows), nCols(nCols), isSorted(true),
          zero(zero) {}

    /**
     * Constructs the @a dimension x @a dimension Matrix from the elements at
     * position @a positions with values @values.
     * @param dimension Defines how many rows and columns this matrix has.
     * @param triplets The nonzero elements.
     * @param zero The zero element (default is 0).
     * @param isSorted True, if the triplets are sorted per row. Default is false.
     */
    CSRGeneralMatrix(count dimension, const std::vector<Triplet> &triplets, ValueType zero = 0,
                     bool isSorted = false)
        : CSRGeneralMatrix(dimension, dimension, triplets, zero, isSorted) {}

    /**
     * Constructs the @a nRows x @a nCols Matrix from the elements at position @a
     * positions with values @values.
     * @param nRows Defines how many rows this matrix has.
     * @param nCols Defines how many columns this matrix has.
     * @param triplets The nonzero elements.
     * @param zero The zero element (default is 0).
     * @param isSorted True, if the triplets are sorted per row. Default is false.
     */
    CSRGeneralMatrix(count nRows, count nCols, const std::vector<Triplet> &triplets,
                     ValueType zero = 0, bool isSorted = false)
        : rowIdx(nRows + 1), columnIdx(triplets.size()), nonZeros(triplets.size()), nRows(nRows),
          nCols(nCols), isSorted(isSorted), zero(zero) {

        const count nnz = triplets.size();

        for (index i = 0; i < nnz; ++i)
            rowIdx[triplets[i].row]++;

        for (index i = 0, prefixSum = 0; i < nRows; ++i) {
            count nnzInRow = rowIdx[i];
            rowIdx[i] = prefixSum;
            prefixSum += nnzInRow;
        }
        rowIdx[nRows] = nnz;

        for (index i = 0; i < nnz; ++i) {
            index row = triplets[i].row;
            index dest = rowIdx[row];

            columnIdx[dest] = triplets[i].column;
            nonZeros[dest] = triplets[i].value;

            rowIdx[row]++;
        }

        for (index i = 0, firstIdxOfRow = 0; i <= nRows; ++i) {
            index newRow = rowIdx[i];
            rowIdx[i] = firstIdxOfRow;
            firstIdxOfRow = newRow;
        }
    }

    /**
     * Constructs the @a nRows x @a nCols Matrix from the elements stored in @a
     * columnIdx and @a values. @a columnIdx and @a values store the colums and
     * values by row.
     * @param nRows
     * @param nCols
     * @param columnIdx
     * @param values
     * @param zero The zero element (default is 0).
     * @param isSorted True if the column indices in @a columnIdx are sorted in
     * every row.
     */
    CSRGeneralMatrix(count nRows, count nCols, const std::vector<std::vector<index>> &columnIdx,
                     const std::vector<std::vector<ValueType>> &values, ValueType zero = 0,
                     bool isSorted = false)
        : rowIdx(nRows + 1), nRows(nRows), nCols(nCols), isSorted(isSorted), zero(zero) {

        count nnz = columnIdx[0].size();
        for (index i = 1; i < columnIdx.size(); ++i) {
            rowIdx[i] = rowIdx[i - 1] + columnIdx[i - 1].size();
            nnz += columnIdx[i].size();
        }
        rowIdx[nRows] = nnz;

        this->columnIdx = std::vector<index>(nnz);
        this->nonZeros = std::vector<double>(nnz);

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
            for (index k = 0; k < columnIdx[i].size(); ++k) {
                this->columnIdx[rowIdx[i] + k] = columnIdx[i][k];
                nonZeros[rowIdx[i] + k] = values[i][k];
            }
        }
    }

    /**
     * Constructs the @a nRows x @a nCols Matrix from the elements at position @a
     * positions with values @values.
     * @param nRows Defines how many rows this matrix has.
     * @param nCols Defines how many columns this matrix has.
     * @param rowIdx The rowIdx vector of the CSR format.
     * @param columnIdx The columnIdx vector of the CSR format.
     * @param nonZeros The nonZero vector of the CSR format. Should be as long as
     * the @a columnIdx vector.
     * @param zero The zero element (default is 0).
     * @param isSorted True, if the triplets are sorted per row. Default is false.
     */
    CSRGeneralMatrix(count nRows, count nCols, const std::vector<index> &rowIdx,
                     const std::vector<index> &columnIdx, const std::vector<ValueType> &nonZeros,
                     ValueType zero = 0, bool isSorted = false)
        : rowIdx(rowIdx), columnIdx(columnIdx), nonZeros(nonZeros), nRows(nRows), nCols(nCols),
          isSorted(isSorted), zero(zero) {}

    /** Default copy constructor */
    CSRGeneralMatrix(const CSRGeneralMatrix &other) = default;

    /** Default move constructor */
    CSRGeneralMatrix(CSRGeneralMatrix &&other) noexcept = default;

    /** Default destructor */
    ~CSRGeneralMatrix() = default;

    /** Default move assignment operator */
    CSRGeneralMatrix &operator=(CSRGeneralMatrix &&other) noexcept = default;

    /** Default copy assignment operator */
    CSRGeneralMatrix &operator=(const CSRGeneralMatrix &other) = default;

    /**
     * Compares this matrix to @a other and returns true if the shape and zero
     * element are the same as well as
     * all entries, otherwise returns false.
     * @param other
     */
    bool operator==(const CSRGeneralMatrix &other) const {
        bool equal = nRows == other.nRows && nCols == other.nCols && zero == other.zero;
        if (equal)
            forNonZeroElementsInRowOrder([&](index i, index j, ValueType value) {
                if (other(i, j) != value) {
                    equal = false;
                    return;
                }
            });

        return equal;
    }

    /**
     * Compares this matrix to @a other and returns false if the shape and zero
     * element are the same as well as
     * all entries, otherwise returns true.
     * @param other
     */
    bool operator!=(const CSRGeneralMatrix &other) const { return !((*this) == other); }

    /**
     * @return Number of rows.
     */
    count numberOfRows() const noexcept { return nRows; }

    /**
     * @return Number of columns.
     */
    count numberOfColumns() const noexcept { return nCols; }

    /**
     * Returns the zero element of the matrix.
     */
    ValueType getZero() const noexcept { return zero; }

    /**
     * @param i The row index.
     * @return Number of non-zeros in row @a i.
     */
    count nnzInRow(const index i) const {
        assert(i < nRows);
        return rowIdx[i + 1] - rowIdx[i];
    }

    /**
     * @return Number of non-zeros in this matrix.
     */
    count nnz() const noexcept { return nonZeros.size(); }

    /**
     * @return Value at matrix position (i,j).
     */
    ValueType operator()(index i, index j) const {
        assert(i < nRows);
        assert(j < nCols);

        if (rowIdx[i] == rowIdx[i + 1])
            return zero; // no non-zero value is present in this row

        double value = zero;
        if (!sorted()) {
            for (index k = rowIdx[i]; k < rowIdx[i + 1]; ++k) {
                if (columnIdx[k] == j) {
                    value = nonZeros[k];
                    break;
                }
            }
        } else {
            index colIdx = binarySearchColumns(rowIdx[i], rowIdx[i + 1] - 1, j);
            if (colIdx != none && rowIdx[i] <= colIdx && columnIdx[colIdx] == j) {
                value = nonZeros[colIdx];
            }
        }

        return value;
    }

    /**
     * Set the matrix at position (@a i, @a j) to @a value.
     * @note This operation can be linear in the number of non-zeros due to vector
     * element movements
     */
    void setValue(index i, index j, ValueType value) {
        assert(i < nRows);
        assert(j < nCols);

        index colIdx = none;
        if (nnzInRow(i) == 0) {
            colIdx = none;
        } else if (!sorted()) {
            for (index k = rowIdx[i]; k < rowIdx[i + 1]; ++k) {
                if (columnIdx[k] == j) {
                    colIdx = k;
                }
            }
        } else {
            colIdx = binarySearchColumns(rowIdx[i], rowIdx[i + 1] - 1, j);
        }

        if (colIdx != none && colIdx >= rowIdx[i]
            && columnIdx[colIdx] == j) { // the matrix already has an entry at (i,j) => replace it
            if (value == getZero()) {    // remove the nonZero value
                columnIdx.erase(columnIdx.begin() + colIdx);
                nonZeros.erase(nonZeros.begin() + colIdx);

                // update rowIdx
                for (index k = i + 1; k < rowIdx.size(); ++k) {
                    --rowIdx[k];
                }
            } else {
                nonZeros[colIdx] = value;
            }
        } else { // create a new non-zero entry at (i,j)
            if (!sorted()) {
                columnIdx.emplace(std::next(columnIdx.begin(), rowIdx[i + 1]), j);
                nonZeros.emplace(std::next(nonZeros.begin(), rowIdx[i + 1]), value);
            } else {
                if (colIdx < rowIdx[i] || colIdx == none) { // emplace the value in
                                                            // front of all other values
                                                            // of row i
                    columnIdx.emplace(std::next(columnIdx.begin(), rowIdx[i]), j);
                    nonZeros.emplace(std::next(nonZeros.begin(), rowIdx[i]), value);
                } else {
                    columnIdx.emplace(std::next(columnIdx.begin(), colIdx + 1), j);
                    nonZeros.emplace(std::next(nonZeros.begin(), colIdx + 1), value);
                }
            }

            // update rowIdx
            for (index k = i + 1; k < rowIdx.size(); ++k) {
                rowIdx[k]++;
            }
        }
    }

    /**
     * Sorts the column indices in each row for faster access.
     */
    void sort() {
#pragma omp parallel for schedule(guided)
        for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
            if (rowIdx[i + 1] - rowIdx[i] > 1) {
                quicksort(rowIdx[i], rowIdx[i + 1] - 1);
            }
        }

#ifndef NETWORKIT_SANITY_CHECKS
        bool sorted = true;
        tlx::unused(sorted);
#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
            for (index j = rowIdx[i] + 1; j < rowIdx[i + 1]; ++j) {
                if (columnIdx[j - 1] > columnIdx[j]) {
                    sorted = false;
                    break;
                }
            }
        }
        assert(sorted);
#endif // NETWORKIT_SANITY_CHECKS

        isSorted = true;
    }

    /**
     * @return True if the matrix is sorted, otherwise false.
     */
    bool sorted() const noexcept { return isSorted; }

    /**
     * @return Row @a i of this matrix as vector.
     */
    Vector row(index i) const {
        assert(i < nRows);

        Vector row(numberOfColumns(), zero, true);
        parallelForNonZeroElementsInRow(i, [&row](index j, double value) { row[j] = value; });

        return row;
    }

    /**
     * @return Column @a j of this matrix as vector.
     */
    Vector column(index j) const {
        assert(j < nCols);

        Vector column(numberOfRows(), getZero());
#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(numberOfRows()); ++i)
            column[i] = (*this)(i, j);

        return column;
    }

    /**
     * @return The main diagonal of this matrix.
     */
    Vector diagonal() const {
        Vector diag(std::min(nRows, nCols), zero);

        if (sorted()) {
#pragma omp parallel for
            for (omp_index i = 0; i < static_cast<omp_index>(diag.getDimension()); ++i) {
                if (rowIdx[i] == rowIdx[i + 1])
                    continue; // no entry in row i
                index left = rowIdx[i];
                index right = rowIdx[i + 1] - 1;
                index mid = (left + right) / 2;
                while (left <= right) {
                    if (columnIdx[mid] == i) {
                        diag[i] = nonZeros[mid];
                        break;
                    }

                    if (columnIdx[mid] < i) {
                        left = mid + 1;
                    } else {
                        right = mid - 1;
                    }

                    mid = (left + right) / 2;
                }
            }
        } else {
#pragma omp parallel for
            for (omp_index i = 0; i < static_cast<omp_index>(diag.getDimension()); ++i) {
                diag[i] = (*this)(i, i);
            }
        }

        return diag;
    }

    /**
     * Adds this matrix to @a other and returns the result.
     * @return The sum of this matrix and @a other.
     */
    CSRGeneralMatrix operator+(const CSRGeneralMatrix &other) const {
        assert(nRows == other.nRows && nCols == other.nCols);
        return CSRGeneralMatrix<ValueType>::binaryOperator(
            *this, other, [](double val1, double val2) { return val1 + val2; });
    }

    /**
     * Adds @a other to this matrix.
     * @return Reference to this matrix.
     */
    CSRGeneralMatrix &operator+=(const CSRGeneralMatrix &other) {
        assert(nRows == other.nRows && nCols == other.nCols);
        *this = CSRGeneralMatrix<ValueType>::binaryOperator(
            *this, other, [](double val1, double val2) { return val1 + val2; });
        return *this;
    }

    /**
     * Subtracts @a other from this matrix and returns the result.
     * @return The difference of this matrix and @a other.
     *
     */
    CSRGeneralMatrix operator-(const CSRGeneralMatrix &other) const {
        assert(nRows == other.nRows && nCols == other.nCols);
        return CSRGeneralMatrix<ValueType>::binaryOperator(
            *this, other, [](double val1, double val2) { return val1 - val2; });
    }

    /**
     * Subtracts @a other from this matrix.
     * @return Reference to this matrix.
     */
    CSRGeneralMatrix &operator-=(const CSRGeneralMatrix &other) {
        assert(nRows == other.nRows && nCols == other.nCols);
        *this = CSRGeneralMatrix<ValueType>::binaryOperator(
            *this, other, [](double val1, double val2) { return val1 - val2; });
        return *this;
    }

    /**
     * Multiplies this matrix with a scalar specified in @a scalar and returns the
     * result.
     * @return The result of multiplying this matrix with @a scalar.
     */
    CSRGeneralMatrix operator*(const ValueType &scalar) const {
        return CSRGeneralMatrix(*this) *= scalar;
    }

    /**
     * Multiplies this matrix with a scalar specified in @a scalar.
     * @return Reference to this matrix.
     */
    CSRGeneralMatrix &operator*=(const ValueType &scalar) {
#pragma omp parallel for
        for (omp_index k = 0; k < static_cast<omp_index>(nonZeros.size()); ++k)
            nonZeros[k] *= scalar;

        return *this;
    }

    /**
     * Multiplies this matrix with @a vector and returns the result.
     * @return The result of multiplying this matrix with @a vector.
     */
    Vector operator*(const Vector &vector) const {
        assert(!vector.isTransposed());
        assert(nCols == vector.getDimension());

        Vector result(nRows, zero);
#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(numberOfRows()); ++i) {
            double sum = zero;
            for (index cIdx = rowIdx[i]; cIdx < rowIdx[i + 1]; ++cIdx) {
                sum += nonZeros[cIdx] * vector[columnIdx[cIdx]];
            }
            result[i] = sum;
        }

        return result;
    }

    /**
     * Multiplies this matrix with @a other and returns the result in a new
     * matrix.
     * @return The result of multiplying this matrix with @a other.
     */
    CSRGeneralMatrix operator*(const CSRGeneralMatrix &other) const {
        assert(nCols == other.nRows);

        std::vector<index> rowIdx(numberOfRows() + 1, 0);
        std::vector<index> columnIdx;
        std::vector<double> nonZeros;

#pragma omp parallel
        {
            std::vector<int64_t> marker(other.numberOfColumns(), -1);
            count numThreads = omp_get_num_threads();
            index threadId = omp_get_thread_num();

            count chunkSize = (numberOfRows() + numThreads - 1) / numThreads;
            index chunkStart = threadId * chunkSize;
            index chunkEnd = std::min(numberOfRows(), chunkStart + chunkSize);

            for (index i = chunkStart; i < chunkEnd; ++i) {
                for (index jA = this->rowIdx[i]; jA < this->rowIdx[i + 1]; ++jA) {
                    index k = this->columnIdx[jA];
                    for (index jB = other.rowIdx[k]; jB < other.rowIdx[k + 1]; ++jB) {
                        index j = other.columnIdx[jB];
                        if (marker[j] != (int64_t)i) {
                            marker[j] = i;
                            ++rowIdx[i + 1];
                        }
                    }
                }
            }

            std::fill(marker.begin(), marker.end(), -1);

#pragma omp barrier
#pragma omp single
            {
                for (index i = 0; i < numberOfRows(); ++i)
                    rowIdx[i + 1] += rowIdx[i];

                columnIdx = std::vector<index>(rowIdx[numberOfRows()]);
                nonZeros = std::vector<double>(rowIdx[numberOfRows()]);
            }

            for (index i = chunkStart; i < chunkEnd; ++i) {
                index rowBegin = rowIdx[i];
                index rowEnd = rowBegin;

                for (index jA = this->rowIdx[i]; jA < this->rowIdx[i + 1]; ++jA) {
                    index k = this->columnIdx[jA];
                    double valA = this->nonZeros[jA];

                    for (index jB = other.rowIdx[k]; jB < other.rowIdx[k + 1]; ++jB) {
                        index j = other.columnIdx[jB];
                        double valB = other.nonZeros[jB];

                        if (marker[j] < (int64_t)rowBegin) {
                            marker[j] = rowEnd;
                            columnIdx[rowEnd] = j;
                            nonZeros[rowEnd] = valA * valB;
                            ++rowEnd;
                        } else {
                            nonZeros[marker[j]] += valA * valB;
                        }
                    }
                }
            }
        }

        CSRGeneralMatrix result(numberOfRows(), other.numberOfColumns(), rowIdx, columnIdx,
                                nonZeros);
        if (sorted() && other.sorted())
            result.sort();
        return result;
    }

    /**
     * Divides this matrix by a divisor specified in @a divisor and returns the
     * result in a new matrix.
     * @return The result of dividing this matrix by @a divisor.
     */
    CSRGeneralMatrix operator/(const ValueType &divisor) const {
        return CSRGeneralMatrix(*this) /= divisor;
    }

    /**
     * Divides this matrix by a divisor specified in @a divisor.
     * @return Reference to this matrix.
     */
    CSRGeneralMatrix &operator/=(const ValueType &divisor) { return *this *= 1.0 / divisor; }

    /**
     * Computes @a A @a binaryOp @a B on the elements of matrix @a A and matrix @a
     * B.
     * @param A Sorted CSRGeneralMatrix.
     * @param B Sorted CSRGeneralMatrix.
     * @param binaryOp Function handling (ValueType, ValueType) -> ValueType
     * @return @a A @a binaryOp @a B.
     * @note @a A and @a B must have the same dimensions and must be sorted.
     */
    template <typename L>
    static CSRGeneralMatrix binaryOperator(const CSRGeneralMatrix &A, const CSRGeneralMatrix &B,
                                           L binaryOp);

    /**
     * Computes @a A^T * @a B.
     * @param A
     * @param B
     * @return @a A^T * @a B.
     * @note The number of rows of @a A must be equal to the number of rows of @a
     * B.
     */
    static CSRGeneralMatrix mTmMultiply(const CSRGeneralMatrix &A, const CSRGeneralMatrix &B) {
        assert(A.nRows == B.nRows);

        std::vector<std::vector<index>> columnIdx(A.numberOfColumns());
        std::vector<std::vector<double>> values(A.numberOfColumns());

        for (index k = 0; k < A.numberOfRows(); ++k) {
            A.forNonZeroElementsInRow(k, [&](index i, double vA) {
                B.forNonZeroElementsInRow(k, [&](index j, double vB) {
                    bool found = false;
                    for (index l = 0; l < columnIdx[i].size(); ++l) {
                        if (columnIdx[i][l] == j) {
                            values[i][l] += vA * vB;
                            found = true;
                            break;
                        }
                    }

                    if (!found) {
                        columnIdx[i].push_back(j);
                        values[i].push_back(vA * vB);
                    }
                });
            });
        }

        return CSRGeneralMatrix(A.nCols, B.nCols, columnIdx, values);
    }

    /**
     * Computes @a A * @a B^T.
     * @param A
     * @param B
     * @return @a A * @a B^T.
     * @note The number of columns of @a A must be equal to the number of columns
     * of @a B.
     */
    static CSRGeneralMatrix mmTMultiply(const CSRGeneralMatrix &A, const CSRGeneralMatrix &B) {
        assert(A.nCols == B.nCols);

        std::vector<std::vector<index>> columnIdx(A.numberOfRows());
        std::vector<std::vector<double>> values(A.numberOfRows());

        for (index i = 0; i < A.numberOfRows(); ++i) {
            A.forNonZeroElementsInRow(i, [&](index k, double vA) {
                for (index j = 0; j < B.numberOfRows(); ++j) {
                    double vB = B(j, k);
                    if (vB != A.zero) {
                        bool found = false;
                        for (index l = 0; l < columnIdx[i].size(); ++l) {
                            if (columnIdx[i][l] == j) {
                                values[i][l] += vA * vB;
                                found = true;
                                break;
                            }
                        }

                        if (!found) {
                            columnIdx[i].push_back(j);
                            values[i].push_back(vA * vB);
                        }
                    }
                }
            });
        }

        return CSRGeneralMatrix(A.nRows, B.nRows, columnIdx, values);
    }

    /**
     * Computes @a matrix^T * @a vector.
     * @param matrix
     * @param vector
     * @return @a matrix^T * @a vector.
     * @note The number of rows of @a matrix must be equal to the dimension of @a
     * vector.
     */
    static Vector mTvMultiply(const CSRGeneralMatrix &matrix, const Vector &vector) {
        assert(matrix.nRows == vector.getDimension() && !vector.isTransposed());

        Vector result(matrix.numberOfColumns(), 0);
        for (index k = 0; k < matrix.numberOfRows(); ++k) {
            matrix.forNonZeroElementsInRow(
                k, [&](index j, double value) { result[j] += value * vector[k]; });
        }

        return result;
    }

    /**
     * Transposes this matrix and returns it.
     */
    CSRGeneralMatrix transpose() const {
        std::vector<index> rowIdx(numberOfColumns() + 1);
        for (index i = 0; i < nnz(); ++i)
            ++rowIdx[columnIdx[i] + 1];

        for (index i = 0; i < numberOfColumns(); ++i)
            rowIdx[i + 1] += rowIdx[i];

        std::vector<index> columnIdx(rowIdx[numberOfColumns()]);
        std::vector<double> nonZeros(rowIdx[numberOfColumns()]);

        for (index i = 0; i < numberOfRows(); ++i) {
            for (index j = this->rowIdx[i]; j < this->rowIdx[i + 1]; ++j) {
                index colIdx = this->columnIdx[j];
                columnIdx[rowIdx[colIdx]] = i;
                nonZeros[rowIdx[colIdx]] = this->nonZeros[j];
                ++rowIdx[colIdx];
            }
        }
        index shift = 0;
        for (index i = 0; i < numberOfColumns(); ++i) {
            index temp = rowIdx[i];
            rowIdx[i] = shift;
            shift = temp;
        }
        rowIdx[numberOfColumns()] = nonZeros.size();

        return CSRGeneralMatrix(nCols, nRows, rowIdx, columnIdx, nonZeros, getZero());
    }

    /**
     * Extracts a matrix with rows and columns specified by @a rowIndices and @a
     * columnIndices from this matrix.
     * The order of rows and columns is equal to the order in @a rowIndices and @a
     * columnIndices. It is also
     * possible to specify a row or column more than once to get duplicates.
     * @param rowIndices
     * @param columnIndices
     */
    CSRGeneralMatrix extract(const std::vector<index> &rowIndices,
                             const std::vector<index> &columnIndices) const {
        std::vector<Triplet> triplets;
        std::vector<std::vector<index>> columnMapping(numberOfColumns());
        for (index j = 0; j < columnIndices.size(); ++j)
            columnMapping[columnIndices[j]].push_back(j);

        bool sorted = true;
        for (index i = 0; i < rowIndices.size(); ++i) {
            Triplet last = {i, 0, 0};
            (*this).forNonZeroElementsInRow(rowIndices[i], [&](index k, double value) {
                if (!columnMapping[k].empty()) {
                    for (index j : columnMapping[k]) {
                        if (last.row == i && last.column > j)
                            sorted = false;
                        last = {i, j, value};
                        triplets.push_back(last);
                    }
                }
            });
        }

        return CSRGeneralMatrix(rowIndices.size(), columnIndices.size(), triplets, getZero(),
                                sorted);
    }

    /**
     * Assign the contents of the matrix @a source to this matrix at rows and
     * columns specified by @a rowIndices and
     * @a columnIndices. That is, entry (i,j) of @a source is assigned to entry
     * (rowIndices[i], columnIndices[j]) of
     * this matrix. Note that the dimensions of @rowIndices and @a columnIndices
     * must coincide with the number of rows
     * and columns of @a source.
     * @param rowIndices
     * @param columnIndices
     * @param source
     */
    void assign(const std::vector<index> &rowIndices, const std::vector<index> &columnIndices,
                const CSRGeneralMatrix &source) {
        assert(rowIndices.size() == source.numberOfRows());
        assert(columnIndices.size() == source.numberOfColumns());

        for (index i = 0; i < rowIndices.size(); ++i)
            source.forElementsInRow(i, [&](index j, double value) {
                setValue(rowIndices[i], columnIndices[j], value);
            });
    }

    /**
     * Applies the unary function @a unaryElementFunction to each value in the
     * matrix. Note that it must hold that the
     * function applied to the zero element of this matrix returns the zero
     * element.
     * @param unaryElementFunction
     */
    template <typename F>
    void apply(F unaryElementFunction);

    /**
     * Compute the (weighted) adjacency matrix of the (weighted) Graph @a graph.
     * @param graph
     */
    static CSRGeneralMatrix adjacencyMatrix(const Graph &graph, ValueType zero = 0) {
        count nonZeros = graph.isDirected() ? graph.numberOfEdges() : graph.numberOfEdges() * 2;
        std::vector<Triplet> triplets(nonZeros);
        index idx = 0;
        graph.forEdges([&](node i, node j, double val) {
            triplets[idx++] = {i, j, val};
            if (!graph.isDirected() && i != j)
                triplets[idx++] = {j, i, val};
        });

        return CSRGeneralMatrix(graph.upperNodeIdBound(), triplets, zero);
    }

    /**
     * Creates a diagonal matrix with dimension equal to the dimension of the
     * Vector @a diagonalElements. The values on
     * the diagonal are the ones stored in @a diagonalElements (i.e. D(i,i) =
     * diagonalElements[i]).
     * @param diagonalElements
     */
    static CSRGeneralMatrix diagonalMatrix(const Vector &diagonalElements, ValueType zero = 0) {
        count nRows = diagonalElements.getDimension();
        count nCols = diagonalElements.getDimension();
        std::vector<index> rowIdx(nRows + 1, 0);
        std::iota(rowIdx.begin(), rowIdx.end(), 0);
        std::vector<index> columnIdx(nCols);
        std::vector<double> nonZeros(nCols);

#pragma omp parallel for
        for (omp_index j = 0; j < static_cast<omp_index>(nCols); ++j) {
            columnIdx[j] = j;
            nonZeros[j] = diagonalElements[j];
        }

        return CSRGeneralMatrix(nRows, nCols, rowIdx, columnIdx, nonZeros, zero);
    }

    /**
     * Returns the (weighted) incidence matrix of the (weighted) Graph @a graph.
     * @param graph
     */
    static CSRGeneralMatrix incidenceMatrix(const Graph &graph, ValueType zero = 0) {
        if (!graph.hasEdgeIds())
            throw std::runtime_error("Graph has no edge Ids. Index edges first by "
                                     "calling graph.indexEdges()");
        std::vector<Triplet> triplets;

        if (graph.isDirected()) {
            graph.forEdges([&](node u, node v, edgeweight weight, edgeid edgeId) {
                if (u != v) {
                    edgeweight w = std::sqrt(weight);
                    triplets.push_back({u, edgeId, w});
                    triplets.push_back({v, edgeId, -w});
                }
            });
        } else {
            graph.forEdges([&](node u, node v, edgeweight weight, edgeid edgeId) {
                if (u != v) {
                    edgeweight w = std::sqrt(weight);
                    if (u < v) { // orientation: small node number -> great node number
                        triplets.push_back({u, edgeId, w});
                        triplets.push_back({v, edgeId, -w});
                    } else {
                        triplets.push_back({u, edgeId, -w});
                        triplets.push_back({v, edgeId, w});
                    }
                }
            });
        }

        return CSRGeneralMatrix(graph.upperNodeIdBound(), graph.upperEdgeIdBound(), triplets, zero);
    }

    /**
     * Compute the (weighted) Laplacian of the (weighted) Graph @a graph.
     * @param graph
     */
    static CSRGeneralMatrix laplacianMatrix(const Graph &graph, ValueType zero = 0) {
        std::vector<Triplet> triples;

        graph.forNodes([&](const index i) {
            double weightedDegree = 0;
            graph.forNeighborsOf(i, [&](const index j, double weight) { // - adjacency matrix
                if (i != j) // exclude diagonal since this would be subtracted by
                            // the adjacency weight
                    weightedDegree += weight;

                triples.push_back({i, j, -weight});
            });

            triples.push_back({i, i, weightedDegree}); // degree matrix
        });

        return CSRGeneralMatrix(graph.upperNodeIdBound(), triples, zero);
    }

    /**
     * Returns the (weighted) normalized Laplacian matrix of the (weighted) Graph
     * @a graph
     * @param graph
     */
    static CSRGeneralMatrix normalizedLaplacianMatrix(const Graph &graph, ValueType zero = 0) {
        std::vector<Triplet> triples;

        std::vector<double> weightedDegrees(graph.upperNodeIdBound(), 0);
        graph.parallelForNodes([&](const node u) { weightedDegrees[u] = graph.weightedDegree(u); });

        graph.forNodes([&](const node i) {
            graph.forNeighborsOf(i, [&](const node j, double weight) {
                if (i != j)
                    triples.push_back(
                        {i, j, -weight / std::sqrt(weightedDegrees[i] * weightedDegrees[j])});
            });

            if (weightedDegrees[i] != 0) {
                if (graph.isWeighted())
                    triples.push_back({i, i, 1 - (graph.weight(i, i)) / weightedDegrees[i]});
                else
                    triples.push_back({i, i, 1});
            }
        });

        return CSRGeneralMatrix(graph.upperNodeIdBound(), triples, zero);
    }

    /**
     * Iterate over all non-zero elements of row @a row in the matrix and call
     * handler(index column, ValueType value)
     */
    template <typename L>
    void forNonZeroElementsInRow(index row, L handle) const {
        for (index k = rowIdx[row]; k < rowIdx[row + 1]; ++k)
            handle(columnIdx[k], nonZeros[k]);
    }

    /**
     * Iterate in parallel over all non-zero elements of row @a row in the matrix
     * and call handler(index column, ValueType value)
     */
    template <typename L>
    void parallelForNonZeroElementsInRow(index row, L handle) const;

    /**
     * Iterate over all elements in row @a i in the matrix and call handle(index
     * column, ValueType value)
     */
    template <typename L>
    void forElementsInRow(index i, L handle) const;

    /**
     * Iterate over all non-zero elements of the matrix in row order and call
     * handler (lambda closure).
     */
    template <typename L>
    void forNonZeroElementsInRowOrder(L handle) const;

    /**
     * Iterate in parallel over all rows and call handler (lambda closure) on
     * non-zero elements of the matrix.
     */
    template <typename L>
    void parallelForNonZeroElementsInRowOrder(L handle) const;
};

template <typename ValueType>
template <typename L>
inline CSRGeneralMatrix<ValueType>
CSRGeneralMatrix<ValueType>::binaryOperator(const CSRGeneralMatrix<ValueType> &A,
                                            const CSRGeneralMatrix<ValueType> &B, L binaryOp) {
    assert(A.nRows == B.nRows && A.nCols == B.nCols);

    if (A.sorted() && B.sorted()) {
        std::vector<index> rowIdx(A.nRows + 1);
        std::vector<std::vector<index>> columns(A.nRows);
        rowIdx[0] = 0;

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(A.nRows); ++i) {
            index k = A.rowIdx[i];
            index l = B.rowIdx[i];
            while (k < A.rowIdx[i + 1] && l < B.rowIdx[i + 1]) {
                if (A.columnIdx[k] < B.columnIdx[l]) {
                    columns[i].push_back(A.columnIdx[k]);
                    ++k;
                } else if (A.columnIdx[k] > B.columnIdx[l]) {
                    columns[i].push_back(B.columnIdx[l]);
                    ++l;
                } else { // A.columnIdx[k] == B.columnIdx[l]
                    columns[i].push_back(A.columnIdx[k]);
                    ++k;
                    ++l;
                }
                ++rowIdx[i + 1];
            }

            while (k < A.rowIdx[i + 1]) {
                columns[i].push_back(A.columnIdx[k]);
                ++k;
                ++rowIdx[i + 1];
            }

            while (l < B.rowIdx[i + 1]) {
                columns[i].push_back(B.columnIdx[l]);
                ++l;
                ++rowIdx[i + 1];
            }
        }

        for (index i = 0; i < A.nRows; ++i)
            rowIdx[i + 1] += rowIdx[i];

        count nnz = rowIdx[A.nRows];
        std::vector<index> columnIdx(nnz);
        std::vector<ValueType> nonZeros(nnz, A.zero);

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(A.nRows); ++i) {
            for (index cIdx = rowIdx[i], j = 0; cIdx < rowIdx[i + 1]; ++cIdx, ++j)
                columnIdx[cIdx] = columns[i][j];

            columns[i].clear();
            columns[i].resize(0);
            columns[i].shrink_to_fit();
        }

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(A.nRows); ++i) {
            index k = A.rowIdx[i];
            index l = B.rowIdx[i];
            for (index cIdx = rowIdx[i]; cIdx < rowIdx[i + 1]; ++cIdx) {
                if (k < A.rowIdx[i + 1] && columnIdx[cIdx] == A.columnIdx[k]) {
                    nonZeros[cIdx] = A.nonZeros[k];
                    ++k;
                }

                if (l < B.rowIdx[i + 1] && columnIdx[cIdx] == B.columnIdx[l]) {
                    nonZeros[cIdx] = binaryOp(nonZeros[cIdx], B.nonZeros[l]);
                    ++l;
                }
            }
        }

        return CSRGeneralMatrix(A.nRows, A.nCols, rowIdx, columnIdx, nonZeros, A.zero, true);
    } else { // A or B not sorted
        std::vector<int64_t> columnPointer(A.nCols, -1);
        std::vector<ValueType> Arow(A.nCols, A.zero);
        std::vector<ValueType> Brow(A.nCols, B.zero);

        std::vector<Triplet> triplets;

        for (index i = 0; i < A.nRows; ++i) {
            index listHead = 0;
            count nnz = 0;

            // search for nonZeros in our own matrix
            for (index k = A.rowIdx[i]; k < A.rowIdx[i + 1]; ++k) {
                index j = A.columnIdx[k];
                Arow[j] = A.nonZeros[k];

                columnPointer[j] = listHead;
                listHead = j;
                nnz++;
            }

            // search for nonZeros in the other matrix
            for (index k = B.rowIdx[i]; k < B.rowIdx[i + 1]; ++k) {
                index j = B.columnIdx[k];
                Brow[j] = B.nonZeros[k];

                if (columnPointer[j]
                    == -1) { // our own matrix does not have a nonZero entry in column j
                    columnPointer[j] = listHead;
                    listHead = j;
                    nnz++;
                }
            }

            // apply operator on the found nonZeros in A and B
            for (count k = 0; k < nnz; ++k) {
                ValueType value = binaryOp(Arow[listHead], Brow[listHead]);
                if (value != A.zero)
                    triplets.push_back({i, listHead, value});

                index temp = listHead;
                listHead = columnPointer[listHead];

                // reset for next row
                columnPointer[temp] = -1;
                Arow[temp] = A.zero;
                Brow[temp] = B.zero;
            }

            nnz = 0;
        }

        return CSRGeneralMatrix(A.numberOfRows(), A.numberOfColumns(), triplets);
    }
}

template <typename ValueType>
template <typename F>
void CSRGeneralMatrix<ValueType>::apply(const F unaryElementFunction) {
#pragma omp parallel for
    for (omp_index k = 0; k < static_cast<omp_index>(nonZeros.size()); ++k)
        nonZeros[k] = unaryElementFunction(nonZeros[k]);
}

template <typename ValueType>
template <typename L>
inline void CSRGeneralMatrix<ValueType>::parallelForNonZeroElementsInRow(index i, L handle) const {
#pragma omp parallel for
    for (omp_index k = rowIdx[i]; k < static_cast<omp_index>(rowIdx[i + 1]); ++k)
        handle(columnIdx[k], nonZeros[k]);
}

template <typename ValueType>
template <typename L>
inline void CSRGeneralMatrix<ValueType>::forElementsInRow(index i, L handle) const {
    Vector rowVector = row(i);
    index j = 0;
    rowVector.forElements([&](ValueType val) { handle(j++, val); });
}

template <typename ValueType>
template <typename L>
inline void CSRGeneralMatrix<ValueType>::forNonZeroElementsInRowOrder(L handle) const {
    for (index i = 0; i < nRows; ++i)
        for (index k = rowIdx[i]; k < rowIdx[i + 1]; ++k)
            handle(i, columnIdx[k], nonZeros[k]);
}

template <typename ValueType>
template <typename L>
inline void CSRGeneralMatrix<ValueType>::parallelForNonZeroElementsInRowOrder(L handle) const {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i)
        for (index k = rowIdx[i]; k < rowIdx[i + 1]; ++k)
            handle(i, columnIdx[k], nonZeros[k]);
}
} /* namespace NetworKit */
#endif // NETWORKIT_ALGEBRAIC_CSR_GENERAL_MATRIX_HPP_
