// no-networkit-format
/*
 * PivotMDS.cpp
 *
 *  Created on: Jul 7, 2016
 *      Author: Michael Wegner
 */

#include <limits>

#include <networkit/auxiliary/PrioQueue.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/viz/PivotMDS.hpp>

namespace NetworKit {

PivotMDS::PivotMDS(const Graph &G, count dim, count numPivots)
    : GraphLayoutAlgorithm(G, dim), dim(dim), numPivots(numPivots) {}

void PivotMDS::run() {
    count n = G->numberOfNodes();
    if (n < numPivots) {
        WARN("Number of Pivots higher than the number of nodes. Setting number of pivots to number of nodes");
        numPivots = n;
    }
    std::vector<node> pivots = GraphTools::randomNodes(*G, numPivots);

    std::vector<Triplet> triplets;
    for (index j = 0; j < numPivots;
         ++j) { // compute distances from pivots to other nodes
        Aux::PrioQueue<edgeweight, node> PQ(n);
        std::vector<edgeweight> dist(n, std::numeric_limits<edgeweight>::max());

        dist[pivots[j]] = 0;
        PQ.insert(0, pivots[j]);
        while (!PQ.empty()) {
            node v = PQ.extractMin().second;
            triplets.push_back({v, j, dist[v]});
            G->forNeighborsOf(v, [&](node w, edgeweight weight) {
                if (dist[v] + weight < dist[w]) {
                    dist[w] = dist[v] + weight;
                    PQ.changeKey(dist[w], w);
                }
            });
        }
    }

    // double center the squared distance matrix
    std::vector<double> rowMean(n, 0.0);
    std::vector<double> colMean(numPivots, 0.0);

    for (Triplet &triplet : triplets) {
        rowMean[triplet.row] += triplet.value / (double)numPivots;
        colMean[triplet.column] += triplet.value / (double)n;
    }

    double grandMean = 0.0;
    double rowDivisor = 1.0 / (2.0 * (double)n);
    double colDivisor = 1.0 / (2.0 * (double)numPivots);
    for (index i = 0; i < n; ++i) {
        grandMean += rowMean[i] * rowDivisor;
    }

    for (index j = 0; j < numPivots; ++j) {
        grandMean += colMean[j] * colDivisor;
    }

    for (Triplet &triplet : triplets) {
        triplet.value = triplet.value - rowMean[triplet.row] -
                        colMean[triplet.column] + grandMean;
    }

    CSRMatrix C(n, numPivots, triplets);

    // compute C^T * C
    CSRMatrix CC = CSRMatrix::mTmMultiply(C, C);
    CC.sort();

    // power iterate to get the first dim largest eigenvectors
    for (index d = 0; d < dim; ++d) {
        Vector eigenvector;
        double eigenvalue;
        powerMethod(CC, numPivots, eigenvector, eigenvalue);

        Vector pos = C * eigenvector;

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
            vertexCoordinates[i][d] = pos[i];
        }

        // remove orthogonal space spanned by the eigenvector
        double sqNorm = eigenvector.length();
        sqNorm *= sqNorm;

        double factor = eigenvalue / sqNorm;

        std::vector<Triplet> triplets;
        for (index i = 0; i < numPivots; ++i) {
            for (index j = 0; j < numPivots; ++j) {
                triplets.push_back(
                    {i, j, factor * eigenvector[i] * eigenvector[j]});
            }
        }

        CSRMatrix eigenMat(numPivots, numPivots, triplets, true);
        CC -= eigenMat;
    }
}

void PivotMDS::powerMethod(const CSRMatrix &mat, count n,
                           Vector &eigenvector, double &eigenvalue) {
    eigenvector = Vector(n);
    for (index i = 0; i < n; ++i) {
        eigenvector[i] = 2.0 * Aux::Random::real() - 1.0;
    }

    Vector old;
    count numIterations = 0;
    do {
        old = eigenvector;
        eigenvector = mat * old;
        eigenvector /= eigenvector.length();

        numIterations++;
    } while ((eigenvector - old).length() > 1e-6 && numIterations < 1500);

    eigenvalue = Vector::innerProduct(mat * eigenvector, eigenvector) /
                 Vector::innerProduct(eigenvector, eigenvector);
}

} /* namespace NetworKit */
