// no-networkit-format
/*
 * MatrixTools.h
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_ALGEBRAIC_MATRIX_TOOLS_HPP_
#define NETWORKIT_ALGEBRAIC_MATRIX_TOOLS_HPP_

#include <atomic>
#include <cmath>

#include <networkit/graph/Graph.hpp>

namespace MatrixTools {

/**
 * Checks if @a matrix is symmetric.
 * @param matrix
 * @return True if @a matrix is symmetric, otherwise false.
 */
template<class Matrix>
bool isSymmetric(const Matrix& matrix) {
    bool output = true;
    matrix.forNonZeroElementsInRowOrder([&] (NetworKit::index i, NetworKit::index j, NetworKit::edgeweight w) {
        if (fabs(matrix(j, i)-w) > NetworKit::FLOAT_EPSILON) {
            output = false;
            return;
        }
    });
    return output;
}

/**
 * Checks if @a matrix is symmetric diagonally dominant (SDD).
 * @param matrix
 * @return True if @a matrix is SDD, false otherwise.
 */
template<class Matrix>
bool isSDD(const Matrix& matrix) {
    if (!isSymmetric(matrix)) {
        return false;
    }

    /* Criterion: a_ii >= \sum_{j != i} a_ij */
    std::vector<double> row_sum(matrix.numberOfRows());
    matrix.parallelForNonZeroElementsInRowOrder([&] (NetworKit::node i, NetworKit::node j, double value) {
        if (i == j) {
            row_sum[i] += value;
        } else {
            row_sum[i] -= fabs(value);
        }
    });

    return std::all_of(row_sum.begin(), row_sum.end(), [] (double val) {return val > -NetworKit::FLOAT_EPSILON;});
}

/**
 * Checks if @a matrix is a Laplacian matrix.
 * @param matrix
 * @return True if @a matrix is a Laplacian matrix, false otherwise.
 */
template<typename Matrix>
bool isLaplacian(const Matrix& matrix) {
    if (!isSymmetric(matrix)) {
        return false;
    }

    /* Criterion: \forall_i \sum_j A_ij = 0  */
    std::vector<double> row_sum(matrix.numberOfRows());
    std::atomic<bool> right_sign(true);
    matrix.parallelForNonZeroElementsInRowOrder([&] (NetworKit::node i, NetworKit::node j, double value) {
        if (i != j && value > NetworKit::FLOAT_EPSILON) {
            right_sign = false;
        }
        row_sum[i] += value;
    });

    return right_sign && std::all_of(row_sum.begin(), row_sum.end(), [] (double val) {return fabs(val) < NetworKit::FLOAT_EPSILON;});
}

/**
 * Computes a graph having the given @a laplacian.
 * @param laplacian
 * @return The graph having a Laplacian equal to @a laplacian.
 */
template<class Matrix>
NetworKit::Graph laplacianToGraph(const Matrix& laplacian) {
    assert(isLaplacian(laplacian));
    NetworKit::Graph G(std::max(laplacian.numberOfRows(), laplacian.numberOfColumns()), true, false);
    laplacian.forNonZeroElementsInRowOrder([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight) {
        if (u != v) { // exclude diagonal
            if (u < v) {
                G.addEdge(u, v, -weight);
            }
        }
    });

    return G;
}

/**
 * Interprets the @a matrix as adjacency matrix of a graph. If @a matrix is non-symmetric, the graph will be directed.
 * @param matrix
 * @return The graph having an adjacency matrix equal to @a matrix.
 */
template<class Matrix>
NetworKit::Graph matrixToGraph(const Matrix& matrix) {
    bool directed = !isSymmetric(matrix);
    NetworKit::Graph G(std::max(matrix.numberOfRows(), matrix.numberOfColumns()), true, directed);
    matrix.forNonZeroElementsInRowOrder([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight) {
        if (directed || u <= v) {
            G.addEdge(u, v, weight);
        }
    });

    return G;
}

}



#endif // NETWORKIT_ALGEBRAIC_MATRIX_TOOLS_HPP_
