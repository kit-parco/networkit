// no-networkit-format
/*
 * DynamicMatrix.cpp
 *
 *  Created on: 13.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <networkit/algebraic/DynamicMatrix.hpp>

namespace NetworKit {

DynamicMatrix::DynamicMatrix() : graph(0, true, true), nRows(0), nCols(0), zero(0.0) {
}

DynamicMatrix::DynamicMatrix(const count dimension, const double zero) : graph(dimension, true, true), nRows(dimension), nCols(dimension), zero(zero) {
}

DynamicMatrix::DynamicMatrix(const count nRows, const count nCols, const double zero) : graph(std::max(nRows, nCols), true, true), nRows(nRows), nCols(nCols), zero(zero) {
}

DynamicMatrix::DynamicMatrix(const count dimension, const std::vector<Triplet>& triplets, const double zero) : DynamicMatrix(dimension, dimension, triplets, zero) {
}

DynamicMatrix::DynamicMatrix(const count nRows, const count nCols, const std::vector<Triplet>& triplets, const double zero) : graph(std::max(nRows, nCols), true, true), nRows(nRows), nCols(nCols), zero(zero) {
    for (size_t k = 0; k < triplets.size(); ++k) {
        assert(triplets[k].row < nRows && triplets[k].column < nCols);
        graph.addEdge(triplets[k].row, triplets[k].column, triplets[k].value);
    }
}

count DynamicMatrix::nnzInRow(const index i) const {
    assert(i < nRows);
    return graph.degree(i);
}

count DynamicMatrix::nnz() const {
    count nnz = 0;
    for (index i = 0; i < nRows; ++i) {
        nnz += nnzInRow(i);
    }

    return nnz;
}

double DynamicMatrix::operator()(const index i, const index j) const {
    assert(i < nRows);
    assert(j < nCols);

    if (!graph.hasEdge(i,j)) return zero;

    return graph.weight(i,j);
}

void DynamicMatrix::setValue(const index i, const index j, const double value) {
    assert(i < nRows);
    assert(j < nCols);

    if (value == getZero() && graph.hasEdge(i,j)) {
        graph.removeEdge(i,j);
    } else {
        graph.setWeight(i, j, value);
    }
}

Vector DynamicMatrix::row(const index i) const {
    assert(i < nRows);

    Vector row(numberOfColumns(), zero, true);
    graph.forEdgesOf(i, [&](node, node j, double value) {
        row[j] = value;
    });

    return row;
}

Vector DynamicMatrix::column(const index j) const {
    assert(j < nCols);

    Vector column(numberOfRows());
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(numberOfRows()); ++i) {
        double val = graph.weight(i,j);
        column[i] = (val == nullWeight? zero : val);
    }

    return column;
}

Vector DynamicMatrix::diagonal() const {
    Vector diag(std::min(nRows, nCols), zero);
    for (index i = 0; i < diag.getDimension(); ++i) {
        diag[i] = (*this)(i,i);
    }

    return diag;
}

DynamicMatrix DynamicMatrix::operator+(const DynamicMatrix &other) const {
    return DynamicMatrix(*this) += other;
}

DynamicMatrix& DynamicMatrix::operator+=(const DynamicMatrix &other) {
    assert(nRows == other.nRows && nCols == other.nCols);

    other.forNonZeroElementsInRowOrder([&](node i, node j, double value) {
        graph.increaseWeight(i, j, value);
    });

    return *this;
}

DynamicMatrix DynamicMatrix::operator-(const DynamicMatrix &other) const {
    return DynamicMatrix(*this) -= other;
}

DynamicMatrix& DynamicMatrix::operator-=(const DynamicMatrix &other) {
    assert(nRows == other.nRows && nCols == other.nCols);

    other.forNonZeroElementsInRowOrder([&](node i, node j, double value) {
        graph.increaseWeight(i, j, -value);
    });

    return *this;
}

DynamicMatrix DynamicMatrix::operator*(const double scalar) const {
    return DynamicMatrix(*this) *= scalar;
}

DynamicMatrix& DynamicMatrix::operator*=(const double scalar) {
    graph.parallelForEdges([&](node i, node j, double value) {
        graph.setWeight(i, j, value * scalar);
    });

    return *this;
}

Vector DynamicMatrix::operator*(const Vector &vector) const {
    assert(!vector.isTransposed());
    assert(nCols == vector.getDimension());
    Vector result(numberOfRows(), zero);

    parallelForNonZeroElementsInRowOrder([&](node i, node j, double value) {
        result[i] += value * vector[j];
    });

    return result;
}

DynamicMatrix DynamicMatrix::operator*(const DynamicMatrix &other) const {
    assert(nCols == other.nRows);

    DynamicMatrix result(numberOfRows(), other.numberOfColumns());
    SparseAccumulator spa(other.numberOfRows());
    for (index r = 0; r < numberOfRows(); ++r) {
        graph.forNeighborsOf(r, [&](node v, double w1){
            other.graph.forNeighborsOf(v, [&](node u, double w2){
                double value = w1 * w2;
                spa.scatter(value, u);
            });
        });

        spa.gather([&](node row, node column, double value){
            result.graph.addEdge(row, column, value);
        });

        spa.increaseRow();
    }

    return result;
}

DynamicMatrix DynamicMatrix::operator/(const double divisor) const {
    return DynamicMatrix(*this) /= divisor;
}

DynamicMatrix& DynamicMatrix::operator/=(const double divisor) {
    return *this *= 1 / divisor;
}

DynamicMatrix DynamicMatrix::mTmMultiply(const DynamicMatrix &A, const DynamicMatrix &B) {
    assert(A.nRows == B.nRows);

    DynamicMatrix C(A.numberOfColumns(), B.numberOfColumns());
    for (index k = 0; k < A.numberOfRows(); ++k) {
        A.graph.forNeighborsOf(k, [&](index i, edgeweight wA) {
            B.graph.forNeighborsOf(k, [&](index j, edgeweight wB) {
                C.graph.increaseWeight(i, j, wA * wB);
            });
        });
    }

    return C;
}

DynamicMatrix DynamicMatrix::mmTMultiply(const DynamicMatrix &A, const DynamicMatrix &B) {
    assert(A.nCols == B.nCols);

    DynamicMatrix C(A.numberOfRows(), B.numberOfRows());
    for (index i = 0; i < A.numberOfRows(); ++i) {
        A.graph.forNeighborsOf(i, [&](index k, edgeweight wA){
            for (index j = 0; j < B.numberOfRows(); ++j) {
                edgeweight wB = B(j,k);
                if (wB != A.zero) {
                    C.graph.increaseWeight(i, j, wA * wB);
                }
            }
        });
    }

    return C;
}

Vector DynamicMatrix::mTvMultiply(const DynamicMatrix &matrix, const Vector &vector) {
    assert(matrix.nRows == vector.getDimension());

    Vector result(matrix.numberOfColumns(), matrix.getZero());
    for (index k = 0; k < matrix.numberOfRows(); ++k) {
        matrix.graph.forNeighborsOf(k, [&](index j, edgeweight w){
            result[j] += w * vector[k];
        });
    }

    return result;
}

DynamicMatrix DynamicMatrix::transpose() const {
    DynamicMatrix transposedMatrix(numberOfColumns(), numberOfRows(), getZero());
    forNonZeroElementsInRowOrder([&](index i, index j, edgeweight weight){
        transposedMatrix.graph.addEdge(j,i,weight);
    });

    return transposedMatrix;
}

DynamicMatrix DynamicMatrix::extract(const std::vector<index>& rowIndices, const std::vector<index>& columnIndices) const {
    std::vector<Triplet> triplets;
    std::vector<std::vector<index>> columnMapping(numberOfColumns());
    for (index j = 0; j < columnIndices.size(); ++j) {
        assert(columnIndices[j] < numberOfColumns());
        columnMapping[columnIndices[j]].push_back(j);
    }

    for (index i = 0; i < rowIndices.size(); ++i) {
        assert(rowIndices[i] < numberOfRows());
        (*this).forNonZeroElementsInRow(rowIndices[i], [&](index k, double value) {
            if (!columnMapping[k].empty()) {
                for (index j : columnMapping[k]) {
                    triplets.push_back({i, j, value});
                }
            }
        });
    }

    return DynamicMatrix(rowIndices.size(), columnIndices.size(), triplets);
}

void DynamicMatrix::assign(const std::vector<index>& rowIndices, const std::vector<index>& columnIndices, const DynamicMatrix& source) {
    assert(rowIndices.size() == source.numberOfRows());
    assert(columnIndices.size() == source.numberOfColumns());

    for (index i = 0; i < rowIndices.size(); ++i) {
        source.forElementsInRow(i, [&](index j, double value) {
            setValue(rowIndices[i], columnIndices[j], value);
        });
    }
}

DynamicMatrix DynamicMatrix::adjacencyMatrix(const Graph& graph, double zero) {
    DynamicMatrix A(graph.upperNodeIdBound(), zero);
    graph.forEdges([&](node u, node v, edgeweight w) {
        A.setValue(u, v, w);
        if (!graph.isDirected()) { // add symmetric value at (v, u)
            A.setValue(v, u, w);
        }
    });

    return A;
}

DynamicMatrix DynamicMatrix::diagonalMatrix(const Vector& diagonalElements, double zero) {
    DynamicMatrix D(diagonalElements.getDimension(), zero);
    for (index i = 0; i < diagonalElements.getDimension(); ++i) {
        D.setValue(i, i, diagonalElements[i]);
    }

    return D;
}

DynamicMatrix DynamicMatrix::incidenceMatrix(const Graph& graph, double zero) {
    if (!graph.hasEdgeIds()) throw std::runtime_error("Graph has no edge Ids. Index edges first by calling graph.indexEdges()");
    DynamicMatrix I(graph.upperNodeIdBound(), graph.upperEdgeIdBound(), zero);
    if (graph.isDirected()) {
        graph.forEdges([&](node u, node v, edgeweight weight, edgeid edgeId) {
            if (u != v) {
                edgeweight w = sqrt(weight);
                I.setValue(u, edgeId, w);
                I.setValue(v, edgeId, -w);
            }
        });
    } else {
        graph.forEdges([&](node u, node v, edgeweight weight, edgeid edgeId){
            if (u != v) {
                edgeweight w = sqrt(weight);
                if (u < v) { // orientation: small node number -> great node number
                    I.setValue(u, edgeId, w);
                    I.setValue(v, edgeId, -w);
                } else {
                    I.setValue(u, edgeId, -w);
                    I.setValue(v, edgeId, w);
                }
            }
        });
    }

    return I;
}

DynamicMatrix DynamicMatrix::laplacianMatrix(const Graph& graph, double zero) {
    DynamicMatrix L(graph.upperNodeIdBound(), zero);
    graph.forNodes([&](const index i){
        double weightedDegree = 0.0;

        graph.forNeighborsOf(i, [&](const index j, double weight) { // - adjacency matrix
            L.setValue(i, j, -weight);
            if (i != j) { // exclude weight of diagonal since it would be subtracted later
                weightedDegree += weight;
            }
        });

        L.setValue(i, i, weightedDegree); // degree matrix
    });

    return L;
}

DynamicMatrix DynamicMatrix::normalizedLaplacianMatrix(const Graph& graph, double zero) {
    DynamicMatrix nL(graph.upperNodeIdBound(), zero);

    std::vector<double> weightedDegrees(graph.upperNodeIdBound(), 0.0);
    graph.parallelForNodes([&](const node u) {
        weightedDegrees[u] = graph.weightedDegree(u);
    });

    graph.forNodes([&](const node i){
        graph.forNeighborsOf(i, [&](const node j, double weight){
            if (i != j) {
                nL.setValue(i, j, -weight/sqrt(weightedDegrees[i] * weightedDegrees[j]));
            }
        });

        if (weightedDegrees[i] != 0.0) {
            if (graph.isWeighted()) {
                nL.setValue(i, i, 1-(graph.weight(i, i)) / weightedDegrees[i]);
            } else {
                nL.setValue(i, i, 1);
            }
        }
    });

    return nL;
}





} /* namespace NetworKit */
