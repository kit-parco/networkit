// no-networkit-format
/*
 * AlgebraicTriangleCounting.h
 *
 *  Created on: Jul 12, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_ALGEBRAIC_ALGORITHMS_ALGEBRAIC_TRIANGLE_COUNTING_HPP_
#define NETWORKIT_ALGEBRAIC_ALGORITHMS_ALGEBRAIC_TRIANGLE_COUNTING_HPP_

#include <networkit/base/Algorithm.hpp>

namespace NetworKit {

/**
 * @ingroup algebraic
 * Implements a triangle counting algorithm for nodes based on algebraic methods.
 */
template<class Matrix>
class AlgebraicTriangleCounting : public Algorithm {
public:
    /**
     * Creates an instance of AlgebraicTriangleCounting for the given Graph @a graph.
     * @param graph
     */
    AlgebraicTriangleCounting(const Graph& graph) : A(Matrix::adjacencyMatrix(graph)), directed(graph.isDirected()) {}

    /**
     * Computes the number of triangles each node is part of. A triangle is considered as a set of nodes (i.e. if there
     * is a triangle (u,v,w) it only counts as one triangle at each node).
     */
    void run() override;

    /**
     * Returns the score of node @a u.
     * @param u
     */
    count score(node u) const {
        assureFinished();
        assert(u < A.numberOfRows());
        return nodeScores[u];
    }

    /**
     * Returns the scores for all nodes of the graph.
     */
    const std::vector<count> &getScores() const {
        assureFinished();
        return nodeScores;
    }

private:
    Matrix A;
    bool directed;
    std::vector<count> nodeScores;
};

template<class Matrix>
void AlgebraicTriangleCounting<Matrix>::run() {
    Matrix powA = A * A * A;

    nodeScores.clear();
    nodeScores.resize(A.numberOfRows(), 0);

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(powA.numberOfRows()); ++i) {
        nodeScores[i] = directed? powA(i,i) : powA(i,i) / 2.0;
    }

    hasRun = true;
}

} /* namespace NetworKit */

#endif // NETWORKIT_ALGEBRAIC_ALGORITHMS_ALGEBRAIC_TRIANGLE_COUNTING_HPP_
