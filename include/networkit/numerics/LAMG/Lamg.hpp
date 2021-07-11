// no-networkit-format
/*
 * Lamg.hpp
 *
 *  Created on: Oct 20, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_NUMERICS_LAMG_LAMG_HPP_
#define NETWORKIT_NUMERICS_LAMG_LAMG_HPP_

#include <vector>
#include <omp.h>

#include <networkit/numerics/LinearSolver.hpp>
#include <networkit/numerics/LAMG/MultiLevelSetup.hpp>
#include <networkit/numerics/LAMG/SolverLamg.hpp>
#include <networkit/numerics/GaussSeidelRelaxation.hpp>
#include <networkit/algebraic/MatrixTools.hpp>
#include <networkit/components/ParallelConnectedComponents.hpp>

namespace NetworKit {

/**
 * @ingroup numerics
 * Represents the interface to the Lean Algebraic Multigrid (LAMG) graph Laplacian linear solver
 * by Oren E. Livne and Achi Brandt.
 * @see Livne, Oren E., and Achi Brandt. "Lean algebraic multigrid (LAMG): Fast graph Laplacian linear solver." SIAM Journal on Scientific Computing 34.4 (2012): B499-B522.
 */
template<class Matrix>
class Lamg : public LinearSolver<Matrix> {
private:
    bool validSetup;
    GaussSeidelRelaxation<Matrix> smoother;
    MultiLevelSetup<Matrix> lamgSetup;
    Matrix laplacianMatrix;
    std::vector<LevelHierarchy<Matrix>> compHierarchies;
    std::vector<SolverLamg<Matrix>> compSolvers;
    std::vector<LAMGSolverStatus> compStati;

    std::vector<Vector> initialVectors;
    std::vector<Vector> rhsVectors;

    count numComponents;
    std::vector<std::vector<index>> components;
    std::vector<index> graph2Components;

    void initializeForOneComponent();

public:
    /**
     * Construct a solver with the given @a tolerance. The relative residual ||Ax-b||/||b|| will be less than or equal to
     * @a tolerance after the solver finished.
     * @param tolerance
     */
    Lamg(const double tolerance = 1e-6) : LinearSolver<Matrix>(tolerance), validSetup(false), lamgSetup(smoother), numComponents(0) {}
    /** Default destructor */
    ~Lamg() override = default;

    /**
     * Compute the multigrid hierarchy for the given Laplacian matrix @a laplacianMatrix.
     * @param laplacianMatrix
     * @note This method also works for disconnected graphs. If you know that the graph is connected,
     * if is faster to use @ref setupConnected instead.
     */
    void setup(const Matrix &laplacianMatrix) override;

    /**
     * Compute the multigrid hierarchy for te given Laplacian matrix @a laplacianMatrix.
     * @param laplacianMatrix
     * @note The graph has to be connected for this method to work. Otherwise the output is undefined.
     */
    void setupConnected(const Matrix &laplacianMatrix) override;

    /**
     * Computes the @a result for the matrix currently setup and the right-hand side @a rhs.
     * The maximum spent time can be specified by @a maxConvergenceTime and the maximum number of iterations can be set
     * by @a maxIterations.
     * @param rhs
     * @param result
     * @param maxConvergenceTime
     * @param maxIterations
     * @return A @ref SolverStatus object which provides some statistics like the final absolute residual.
     */
    SolverStatus solve(const Vector &rhs, Vector &result, count maxConvergenceTime = 5 * 60 * 1000,
                       count maxIterations = std::numeric_limits<count>::max()) override;

    /**
     * Compute the @a results for the matrix currently setup and the right-hand sides @a rhs.
     * The maximum spent time for each system can be specified by @a maxConvergenceTime and the maximum number of iterations can be set
     * by @a maxIterations.
     * @param rhs
     * @param results
     * @param maxConvergenceTime
     * @param maxIterations
     */
    void parallelSolve(const std::vector<Vector> &rhs, std::vector<Vector> &results,
                       count maxConvergenceTime = 5 * 60 * 1000,
                       count maxIterations = std::numeric_limits<count>::max()) override;

    /**
     * Abstract parallel solve function that computes and processes results using @a resultProcessor for the matrix currently setup and the right-hand sides (size of @a rhsSize) provided by @a rhsLoader.
     * The maximum spent time for each system can be specified by @a maxConvergenceTime and the maximum number of iterations can be set
     * by @a maxIterations.
     * @param rhsLoader
     * @param resultProcessor
     * @param rhsSize
     * @param maxConvergenceTime
     * @param maxIterations
     * @note If the solver does not support parallelism during solves, this function falls back to solving the systems sequentially.
     */
    template<typename RHSLoader, typename ResultProcessor>
    void parallelSolve(const RHSLoader& rhsLoader, const ResultProcessor& resultProcessor, std::pair<count, count> rhsSize,
                       count maxConvergenceTime = 5 * 60 * 1000,
                       count maxIterations = std::numeric_limits<count>::max()) {
        if (numComponents == 1) {
            const index numThreads = omp_get_max_threads();
            if (compSolvers.size() != numThreads) {
                compSolvers.clear();

                for (index i = 0; i < (index) numThreads; ++i) {
                    compSolvers.push_back(SolverLamg<Matrix>(compHierarchies[0], smoother));
                }
            }

            count n = rhsSize.first;
            count m = rhsSize.second;
            std::vector<Vector> results(numThreads, Vector(m));
            std::vector<Vector> RHSs(numThreads, Vector(m));

#pragma omp parallel for
            for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
                index threadId = omp_get_thread_num();

                const Vector& rhs = rhsLoader(i, RHSs[threadId]);
                Vector& result = results[threadId];

                LAMGSolverStatus stat;
                stat.desiredResidualReduction = this->tolerance * rhs.length() / (laplacianMatrix * result - rhs).length();
                stat.maxIters = maxIterations;
                stat.maxConvergenceTime = maxConvergenceTime;

                compSolvers[threadId].solve(result, rhs, stat);
                resultProcessor(i, result);
            }
        }
    }
};

template<class Matrix>
void Lamg<Matrix>::initializeForOneComponent() {
    compHierarchies = std::vector<LevelHierarchy<Matrix>>(1);
    lamgSetup.setup(laplacianMatrix, compHierarchies[0]);
    compSolvers.clear();
    compSolvers.push_back(SolverLamg<Matrix>(compHierarchies[0], smoother));
    validSetup = true;
}

template<class Matrix>
void Lamg<Matrix>::setupConnected(const Matrix& laplacianMatrix) {
    this->laplacianMatrix = laplacianMatrix;
    initializeForOneComponent();
    numComponents = 1;
}

template<class Matrix>
void Lamg<Matrix>::setup(const Matrix& laplacianMatrix) {
    this->laplacianMatrix = laplacianMatrix;
    Graph G = MatrixTools::matrixToGraph(laplacianMatrix);
    ParallelConnectedComponents con(G, false);
    con.run();
    numComponents = con.numberOfComponents();
    if (numComponents == 1) {
        initializeForOneComponent();
    } else {
        graph2Components = std::vector<index>(G.numberOfNodes());

        initialVectors = std::vector<Vector>(numComponents);
        rhsVectors = std::vector<Vector>(numComponents);

        components = std::vector<std::vector<index>>(numComponents);
        compHierarchies = std::vector<LevelHierarchy<Matrix>>(numComponents);
        compSolvers.clear();
        compStati = std::vector<LAMGSolverStatus>(numComponents);

        // create solver for every component
        index compIdx = 0;
        for (const auto &component : con.getPartition().getSubsets()) {
            components[compIdx] = std::vector<index>(component.begin(), component.end());

            std::vector<Triplet> triplets;
            index idx = 0;
            for (node u : components[compIdx]) {
                graph2Components[u] = idx;
                idx++;
            }

            for (node u : components[compIdx]) {
                G.forNeighborsOf(u, [&](node v, edgeweight w) {
                    triplets.push_back({graph2Components[u], graph2Components[v], w});
                });
            }

            Matrix compMatrix(component.size(), component.size(), triplets);
            initialVectors[compIdx] = Vector(component.size());
            rhsVectors[compIdx] = Vector(component.size());
            lamgSetup.setup(compMatrix, compHierarchies[compIdx]);
            compSolvers.push_back(SolverLamg<Matrix>(compHierarchies[compIdx], smoother));
            LAMGSolverStatus status{};
            status.desiredResidualReduction = this->tolerance * component.size() / G.numberOfNodes();
            compStati[compIdx] = status;

            compIdx++;
        }

        validSetup = true;
    }
}

template<class Matrix>
SolverStatus Lamg<Matrix>::solve(const Vector& rhs, Vector& result, count maxConvergenceTime, count maxIterations) {
    if (!validSetup || result.getDimension() != laplacianMatrix.numberOfColumns()
            || rhs.getDimension() != laplacianMatrix.numberOfRows()) {
        throw std::runtime_error("No or wrong matrix is setup for given vectors.");
    }

    SolverStatus status;

    if (numComponents == 1) {
        LAMGSolverStatus stat;
        stat.desiredResidualReduction = this->tolerance * rhs.length() / (laplacianMatrix * result - rhs).length();
        stat.maxIters = maxIterations;
        stat.maxConvergenceTime = maxConvergenceTime;
        compSolvers[0].solve(result, rhs, stat);

        status.residual = stat.residual;
        status.numIters = stat.numIters;
        status.converged = stat.converged;
    } else {
        // solve on every component
        count maxIters = 0;
        for (index i = 0; i < components.size(); ++i) {
            for (auto element : components[i]) {
                initialVectors[i][graph2Components[element]] = result[element];
                rhsVectors[i][graph2Components[element]] = rhs[element];
            }

            double resReduction = this->tolerance * rhsVectors[i].length() / (compHierarchies[i].at(0).getLaplacian() * initialVectors[i] - rhsVectors[i]).length();
            compStati[i].desiredResidualReduction = resReduction * components[i].size() / laplacianMatrix.numberOfRows();
            compStati[i].maxIters = maxIterations;
            compStati[i].maxConvergenceTime = maxConvergenceTime;
            compSolvers[i].solve(initialVectors[i], rhsVectors[i], compStati[i]);

            for (auto element : components[i]) { // write solution back to result
                result[element] = initialVectors[i][graph2Components[element]];
            }

            maxIters = std::max(maxIters, compStati[i].numIters);
        }

        status.residual = (rhs - laplacianMatrix * result).length();
        status.converged = status.residual <= this->tolerance;
        status.numIters = maxIters;
    }

    return status;
}

template<class Matrix>
void Lamg<Matrix>::parallelSolve(const std::vector<Vector>& rhs, std::vector<Vector>& results, count maxConvergenceTime, count maxIterations) {
    if (numComponents == 1) {
        assert(rhs.size() == results.size());
        const index numThreads = omp_get_max_threads();
        if (compSolvers.size() != numThreads) {
            compSolvers.clear();

            for (index i = 0; i < (index) numThreads; ++i) {
                compSolvers.push_back(SolverLamg<Matrix>(compHierarchies[0], smoother));
            }
        }

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(rhs.size()); ++i) {
            index threadId = omp_get_thread_num();
            LAMGSolverStatus stat;
            stat.desiredResidualReduction = this->tolerance * rhs[i].length() / (laplacianMatrix * results[i] - rhs[i]).length();
            stat.maxIters = maxIterations;
            stat.maxConvergenceTime = maxConvergenceTime;
            compSolvers[threadId].solve(results[i], rhs[i], stat);
        }

    }
}

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_LAMG_LAMG_HPP_
