// no-networkit-format
/*
 * AlgebraicPageRank.h
 *
 *  Created on: Jun 20, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_ALGEBRAIC_ALGORITHMS_ALGEBRAIC_PAGE_RANK_HPP_
#define NETWORKIT_ALGEBRAIC_ALGORITHMS_ALGEBRAIC_PAGE_RANK_HPP_

#include <networkit/centrality/Centrality.hpp>
#include <networkit/auxiliary/Parallel.hpp>

#include <networkit/graph/Graph.hpp>
#include <networkit/algebraic/GraphBLAS.hpp>
#include <networkit/algebraic/Vector.hpp>

namespace NetworKit {

/**
 * @ingroup algebraic
 * Implementation of PageRank using the GraphBLAS interface.
 */
template<class Matrix>
class AlgebraicPageRank final : public Centrality {
public:

    /**
     * Constructs an instance of AlgebraicPageRank for the given @a graph. Page rank uses the damping factor @a damp
     * and the tolerance @a tol.
     * @param graph
     * @param damp
     * @param tol
     */
    AlgebraicPageRank(const Graph& graph, const double damp = 0.85, const double tol = 1e-8)
        : Centrality(graph), damp(damp), tol(tol) {
        Matrix A = Matrix::adjacencyMatrix(graph);
        // normalize At by out-degree
        Vector invOutDeg = GraphBLAS::rowReduce(A);
#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(invOutDeg.getDimension()); ++i) {
            invOutDeg[i] = 1.0/invOutDeg[i];
        }

        std::vector<Triplet> mTriplets(A.nnz());
        index idx = 0;
        A.forNonZeroElementsInRowOrder([&](index i, index j, double value) {
            mTriplets[idx++] = {j,i, damp * value * invOutDeg[i]};
        });
        M = std::move(Matrix(A.numberOfRows(), mTriplets));
    }

    void run() override;

    /**
     * Get a vector containing the betweenness score for each node in the graph.
     *
     * @return The betweenness scores calculated by @ref run().
     */
    const std::vector<double> &scores() const override;

    /**
     * Get a vector of pairs sorted into descending order. Each pair contains a node and the corresponding score
     * calculated by @ref run().
     * @return A vector of pairs.
     */
    std::vector<std::pair<node, double>> ranking() override;

    /**
     * Get the betweenness score of node @a v calculated by @ref run().
     *
     * @param v A node.
     * @return The betweenness score of node @a v.
     */
    double score(node v) override;

    /**
     * Get the theoretical maximum of centrality score in the given graph.
     *
     * @return The maximum centrality score.
     */
    double maximum() override {
        return 1.0;
    }

private:
    Matrix M;

    const double damp;
    const double tol;
};

template<class Matrix>
void AlgebraicPageRank<Matrix>::run() {
    count n = M.numberOfRows();
    double teleportProb = (1.0 - damp) / (double) n;
    Vector rank(n, 1.0/(double)n);
    Vector lastRank;

    do {
        lastRank = rank;
        rank = M * rank;
        rank.apply([&](double value) {return value += teleportProb;});
    } while ((rank - lastRank).length() > tol);

    double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for (omp_index i = 0; i < static_cast<omp_index>(rank.getDimension()); ++i) {
        sum += rank[i];
    }

    scoreData.resize(n, 0);
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(rank.getDimension()); ++i) {
        scoreData[i] = rank[i] / sum;
    }

    hasRun = true;
}

template<class Matrix>
const std::vector<double> &AlgebraicPageRank<Matrix>::scores() const {
    assureFinished();
    return scoreData;
}

template<class Matrix>
std::vector<std::pair<node, double>> AlgebraicPageRank<Matrix>::ranking() {
    assureFinished();
    std::vector<std::pair<node, double> > ranking;
    for (index i = 0; i < scoreData.size(); ++i) {
        ranking.push_back({i, scoreData[i]});
    }
    Aux::Parallel::sort(ranking.begin(), ranking.end(), [](std::pair<node, double> x, std::pair<node, double> y) { return x.second > y.second; });
    return ranking;
}

template<class Matrix>
double AlgebraicPageRank<Matrix>::score(node v) {
    assureFinished();
    return scoreData.at(v);
}

} /* namespace NetworKit */

#endif // NETWORKIT_ALGEBRAIC_ALGORITHMS_ALGEBRAIC_PAGE_RANK_HPP_
