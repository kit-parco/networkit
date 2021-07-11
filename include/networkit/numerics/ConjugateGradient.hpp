// no-networkit-format
/*
 * ConjugateGradient.h
 *
 *  Created on: 15.06.2014
 *      Author: Daniel Hoske and Michael Wegner
 */

#ifndef NETWORKIT_NUMERICS_CONJUGATE_GRADIENT_HPP_
#define NETWORKIT_NUMERICS_CONJUGATE_GRADIENT_HPP_

#include <cstdint>
#include <utility>

#include <networkit/numerics/LinearSolver.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>

namespace NetworKit {

/**
 * @ingroup numerics
 * Implementation of Conjugate Gradient.
 */
template<class Matrix, class Preconditioner>
class ConjugateGradient : public LinearSolver<Matrix> {
public:
    ConjugateGradient(double tolerance = 1e-5) : LinearSolver<Matrix>(tolerance), matrix(Matrix()) {}

    void setup(const Matrix &matrix) override {
        this->matrix = matrix;
        precond = Preconditioner(matrix);
    }

    void setupConnected(const Matrix &matrix) override {
        this->matrix = matrix;
        precond = Preconditioner(matrix);
    }

    /**
     * Solves the linear system \f$Ax = b\f$ using the conjugate gradient method
     * with a given preconditioner and with initial value \f$(0, \dots, 0)^T\f$.
     * We the return the solution \f$x\f$. The solution \f$x\f$ fulfils
     * \f$\frac{\Vert Ax - b\Vert}{\Vert b \Vert} \leq relative\_residual\f$ if the
     * algorithm has converged.
     *
     * Obviously, @a A needs to have the same number of rows as @a b and
     * @a status.residual must be nonnegative. You may also request that the algorithm
     * does not run for more than @a status.max_iters iterations.
     */
    SolverStatus solve(const Vector &rhs, Vector &result, count maxConvergenceTime = 5 * 60 * 1000,
                       count maxIterations = std::numeric_limits<count>::max()) override;

    /**
     * Solves the linear systems in parallel.
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
        const index numThreads = omp_get_max_threads();
        count n = rhsSize.first;
        count m = rhsSize.second;
        std::vector<Vector> results(numThreads, Vector(m));
        std::vector<Vector> RHSs(numThreads, Vector(m));

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
            index threadId = omp_get_thread_num();
            const Vector& rhs = rhsLoader(i, RHSs[threadId]);
            Vector& result = results[threadId];

            this->solve(rhs, result, maxConvergenceTime, maxIterations);
            resultProcessor(i, result);
        }
    }

private:
    Matrix matrix;
    Preconditioner precond;

};

template<class Matrix, class Preconditioner>
SolverStatus ConjugateGradient<Matrix, Preconditioner>::solve(const Vector& rhs, Vector& result, count, count maxIterations) {
    assert(matrix.numberOfRows() == rhs.getDimension());

    // Absolute residual to achieve
    double sqr_desired_residual = this->tolerance * this->tolerance * (rhs.length() * rhs.length());

    // Main loop. See: http://en.wikipedia.org/wiki/Conjugate_gradient_method#The_resulting_algorithm
    Vector residual_dir = rhs - matrix*result;
    Vector conjugate_dir = precond.rhs(residual_dir);
    double sqr_residual = Vector::innerProduct(residual_dir, residual_dir);
    double sqr_residual_precond = Vector::innerProduct(residual_dir, conjugate_dir);

    count niters = 0;
    Vector tmp, residual_precond;
    while (sqr_residual > sqr_desired_residual) {
        niters++;
        if (niters > maxIterations) {
            break;
        }

        tmp = matrix * conjugate_dir;
        double step = sqr_residual_precond / Vector::innerProduct(conjugate_dir, tmp);
        result += step * conjugate_dir;
        residual_dir -= step * tmp;
        sqr_residual = Vector::innerProduct(residual_dir, residual_dir);

        residual_precond = precond.rhs(residual_dir);
        double new_sqr_residual_precond = Vector::innerProduct(residual_dir, residual_precond);
        conjugate_dir = (new_sqr_residual_precond / sqr_residual_precond) * conjugate_dir + residual_precond;
        sqr_residual_precond = new_sqr_residual_precond;
    }

    SolverStatus status;
    status.numIters = niters;
    status.residual = (rhs - matrix*result).length();
    status.converged = status.residual / rhs.length() <= this->tolerance;

    return status;
}

template<class Matrix, class Preconditioner>
void ConjugateGradient<Matrix, Preconditioner>::parallelSolve(const std::vector<Vector>& rhs, std::vector<Vector>& results, count maxConvergenceTime, count maxIterations) {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(rhs.size()); ++i) {
        this->solve(rhs[i], results[i], maxConvergenceTime, maxIterations);
    }
}

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_CONJUGATE_GRADIENT_HPP_
