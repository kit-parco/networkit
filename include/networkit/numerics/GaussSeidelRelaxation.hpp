// no-networkit-format
/*
 * GaussSeidelRelaxation.h
 *
 *  Created on: 27.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_NUMERICS_GAUSS_SEIDEL_RELAXATION_HPP_
#define NETWORKIT_NUMERICS_GAUSS_SEIDEL_RELAXATION_HPP_

#include <networkit/numerics/Smoother.hpp>

namespace NetworKit {

/**
 * @ingroup numerics
 * Implementation of the Gauss-Seidel smoother.
 */
template<class Matrix>
class GaussSeidelRelaxation : public Smoother<Matrix> {

private:
    double tolerance;

public:
    /**
     * Constructs a Gauss-Seidel smoother with the given @a tolerance (default: 1e-15).
     * @param tolerance
     */
    GaussSeidelRelaxation(double tolerance=1e-15) : tolerance(tolerance) {}

    /**
     * Utilizes Gauss-Seidel relaxations until the given number of @a maxIterations is reached or the relative residual
     * is below the tolerance specified in the constructor. The solver starts with @a initialGuess as intitial guess to
     * the solution.
     * @param A The matrix.
     * @param b The right-hand-side.
     * @param initialGuess
     * @param maxIterations
     * @return The (approximate) solution to the system.
     */
    Vector relax(const Matrix &A, const Vector &b, const Vector &initialGuess,
                 count maxIterations = std::numeric_limits<count>::max()) const override;

    /**
     * Utilizes Gauss-Seidel relaxations until the given number of @a maxIterations is reached or the relative residual
     * is below the tolerance specified in the constructor.
     * @param A The matrix.
     * @param b The right-hand-side.
     * @param maxIterations
     * @return The (approximate) solution to the system.
     */
    Vector relax(const Matrix &A, const Vector &b,
                 count maxIterations = std::numeric_limits<count>::max()) const override;
};

template<class Matrix>
Vector GaussSeidelRelaxation<Matrix>::relax(const Matrix& A, const Vector& b, const Vector& initialGuess, const count maxIterations) const {
    count iterations = 0;
    Vector x_old = initialGuess;
    Vector x_new = initialGuess;
    if (maxIterations == 0) return initialGuess;

    count dimension = A.numberOfColumns();
    Vector diagonal = A.diagonal();

    do {
        x_old = x_new;

        for (index i = 0; i < dimension; ++i) {
            double sigma = 0.0;
            A.forNonZeroElementsInRow(i, [&](index column, double value) {
                if (column != i) {
                    sigma += value * x_new[column];
                }
            });

            x_new[i] = (b[i] - sigma) / diagonal[i];
        }

        iterations++;
    } while (iterations < maxIterations && (A*x_new - b).length() / b.length() > tolerance);

    return x_new;
}

template<class Matrix>
Vector GaussSeidelRelaxation<Matrix>::relax(const Matrix& A, const Vector& b, const count maxIterations) const {
    Vector x(b.getDimension());
    return relax(A, b, x, maxIterations);
}


} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_GAUSS_SEIDEL_RELAXATION_HPP_
