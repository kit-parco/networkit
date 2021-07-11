// no-networkit-format
/*
 * SolverLamg.h
 *
 *  Created on: 12.01.2015
 *      Author: Michael
 */

#ifndef NETWORKIT_NUMERICS_LAMG_SOLVER_LAMG_HPP_
#define NETWORKIT_NUMERICS_LAMG_SOLVER_LAMG_HPP_

#include <cmath>
#include <vector>

#include <networkit/algebraic/DenseMatrix.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/numerics/LAMG/LevelHierarchy.hpp>
#include <networkit/numerics/Smoother.hpp>

namespace NetworKit {

/**
 * Status parameters of the solver.
 */
struct LAMGSolverStatus {
    // in
    count maxIters = std::numeric_limits<count>::max(); // maximum number of iterations
    count maxConvergenceTime = std::numeric_limits<count>::max(); // maximum time in milliseconds spent to solve the system
    double desiredResidualReduction = 1e-8; // desired reduction of the initial residual (finalResidual <= desiredResReduction * initialResidual)
    count numPreSmoothIters = 1; // number of pre smoothing iterations
    count numPostSmoothIters = 2; // number of post smoothing iterations

    // out
    count numIters; // number of iterations needed during solve phase
    double residual; // absolute final residual
    bool converged; // flag of conversion status
    std::vector<double> residualHistory; // history of absolute residuals
};

/**
 * @ingroup numerics
 * Implements the solve phase of LAMG (Lean Algebraic Multigrid by Livne et al.).
 */
template<class Matrix>
class SolverLamg {
private:
    LevelHierarchy<Matrix> &hierarchy;
    const Smoother<Matrix> &smoother;

    // data structures for iterate recombination
    std::vector<std::vector<Vector>> history;
    std::vector<std::vector<Vector>> rHistory;
    std::vector<index> latestIterate;
    std::vector<count> numActiveIterates;

    // bStages for Elimination Levels
    std::vector<std::vector<Vector>> bStages;

    void solveCycle(Vector& x, const Vector& b, int finest, LAMGSolverStatus& status);
    void cycle(Vector& x, const Vector& b, int finest, int coarsest, std::vector<count>& numVisits, std::vector<Vector>& X, std::vector<Vector>& B, const LAMGSolverStatus& status);
    void multigridCycle(index level, Vector& xf, const Vector& bf);
    void saveIterate(index level, const Vector& x, const Vector& r);
    void clearHistory(index level);
    void minRes(index level, Vector& x, const Vector& r);

public:
    /**
     * Constructs a new solver instance for the specified @a hierarchy. The @a smoother will be used for relaxing and
     * solving the coarser solutions.
     * @param hierarchy Reference to the LevelHierarchy constructed by MultiLevelSetup.
     * @param smoother Reference to a smoother.
     */
    SolverLamg(LevelHierarchy<Matrix>& hierarchy, const Smoother<Matrix>& smoother) : hierarchy(hierarchy), smoother(smoother), bStages(hierarchy.size(), std::vector<Vector>()) {}

    SolverLamg (const SolverLamg<Matrix>& other) = default;

    SolverLamg (SolverLamg<Matrix>&& other) noexcept = default;

    virtual ~SolverLamg() = default;

    SolverLamg& operator=(SolverLamg<Matrix>&& other) noexcept = default;

    SolverLamg& operator=(const SolverLamg<Matrix>& other) = default;

    /**
     * Solves the system A*x = b for the given initial @a x and right-hand side @a b. More parameters can be specified
     * in @a status and additional output is also stored in @a status. After the solver finished, the approximate
     * solution is stored in @a x.
     * @param x[out] Reference to the initial guess to the solution and the approximation after the solver finished.
     * @param b The right-hand side vector.
     * @param status Reference to an LAMGSolverStatus.
     */
    void solve(Vector& x, const Vector& b, LAMGSolverStatus& status);
};

template<class Matrix>
void SolverLamg<Matrix>::solve(Vector& x, const Vector& b, LAMGSolverStatus& status) {
    bStages = std::vector<std::vector<Vector>>(hierarchy.size(), std::vector<Vector>());
    if (hierarchy.size() >= 2) {
        Vector bc = b;
        Vector xc = x;
        int finest = 0;

        if (hierarchy.getType(1) == ELIMINATION) {
            hierarchy.at(1).restrict(b, bc, bStages[1]);
            if (hierarchy.at(1).getLaplacian().numberOfRows() == 1) {
                x = 0.0;
                return;
            } else {
                hierarchy.at(1).coarseType(x, xc);
                finest = 1;
            }
        }
        solveCycle(xc, bc, finest, status);

        if (finest == 1) { // interpolate from finest == ELIMINATION level back to actual finest level
            hierarchy.at(1).interpolate(xc, x, bStages[1]);
        } else {
            x = xc;
        }
    } else {
        solveCycle(x, b, 0, status);
    }

    double residual = (b - hierarchy.at(0).getLaplacian() * x).length();
    status.residual = residual;
}

template<class Matrix>
void SolverLamg<Matrix>::solveCycle(Vector& x, const Vector& b, int finest, LAMGSolverStatus& status) {
    Aux::Timer timer;
    timer.start();

    // data structures for iterate recombination
    history = std::vector<std::vector<Vector>>(hierarchy.size());
    rHistory = std::vector<std::vector<Vector>>(hierarchy.size());
    latestIterate = std::vector<index>(hierarchy.size(), 0);
    numActiveIterates = std::vector<count>(hierarchy.size(), 0);
    int coarsest = hierarchy.size() - 1;
    std::vector<count> numVisits(coarsest);
    std::vector<Vector> X(hierarchy.size());
    std::vector<Vector> B(hierarchy.size());

    for (index i = 0; i < hierarchy.size(); ++i) {
        history[i] = std::vector<Vector>(MAX_COMBINED_ITERATES, Vector(hierarchy.at(i).getNumberOfNodes()));
        rHistory[i] = std::vector<Vector>(MAX_COMBINED_ITERATES, Vector(hierarchy.at(i).getNumberOfNodes()));
    }

    Vector r = b - hierarchy.at(finest).getLaplacian() * x;
    double residual = r.length();
    double finalResidual = residual * status.desiredResidualReduction;
    double bestResidual = std::numeric_limits<double>::max();

    count iterations = 0;
    status.residualHistory.emplace_back(residual);
    count noResReduction = 0;
    while (residual > finalResidual && noResReduction < 5 && iterations < status.maxIters && timer.elapsedMilliseconds() <= status.maxConvergenceTime ) {
        cycle(x, b, finest, coarsest, numVisits, X, B, status);
        r = b - hierarchy.at(finest).getLaplacian() * x;
        residual = r.length();
        status.residualHistory.emplace_back(residual);
        if (residual < bestResidual) {
            noResReduction = 0;
            bestResidual = residual;
        } else {
            ++noResReduction;
        }
        iterations++;
    }

    timer.stop();

    status.numIters = iterations;
    status.residual = r.length();
    status.converged = r.length() <= finalResidual;
}

template<class Matrix>
void SolverLamg<Matrix>::cycle(Vector& x, const Vector& b, int finest, int coarsest, std::vector<count>& numVisits, std::vector<Vector>& X, std::vector<Vector>& B, const LAMGSolverStatus& status) {
    std::fill(numVisits.begin(), numVisits.end(), 0);
    X[finest] = x;
    B[finest] = b;

    int currLvl = finest;
    int nextLvl = finest;
    double maxVisits = 0.0;

    saveIterate(currLvl, X[currLvl], B[currLvl] - hierarchy.at(currLvl).getLaplacian() * X[currLvl]);
    while (true) {
        if (currLvl == coarsest) {
            nextLvl = currLvl - 1;
            if (currLvl == finest) { // finest level
                X[currLvl] = smoother.relax(hierarchy.at(currLvl).getLaplacian(), B[currLvl], X[currLvl], status.numPreSmoothIters);
            } else {
                Vector bCoarse(B[currLvl].getDimension()+1, 0.0);
                for (index i = 0; i < B[currLvl].getDimension(); ++i) {
                    bCoarse[i] = B[currLvl][i];
                }

                Vector xCoarse = DenseMatrix::LUSolve(hierarchy.getCoarseMatrix(), bCoarse);
                for (index i = 0; i < X[currLvl].getDimension(); ++i) {
                    X[currLvl][i] = xCoarse[i];
                }
            }
        } else {
            if (currLvl == finest) {
                maxVisits = 1.0;
            } else {
                maxVisits = hierarchy.cycleIndex(currLvl) * numVisits[currLvl-1];
            }

            if (numVisits[currLvl] < static_cast<count>(maxVisits)) {
                nextLvl = currLvl + 1;
            } else {
                nextLvl = currLvl - 1;
            }
        }

        if (nextLvl < finest) break;

        if (nextLvl > currLvl) {  // preProcess
            numVisits[currLvl]++;

            if (hierarchy.getType(nextLvl) != ELIMINATION) {
                X[currLvl] = smoother.relax(hierarchy.at(currLvl).getLaplacian(), B[currLvl], X[currLvl], status.numPreSmoothIters);
            }

            if (hierarchy.getType(nextLvl) == ELIMINATION) {
                hierarchy.at(nextLvl).restrict(B[currLvl], B[nextLvl], bStages[nextLvl]);
            } else {
                hierarchy.at(nextLvl).restrict(B[currLvl] - hierarchy.at(currLvl).getLaplacian() * X[currLvl], B[nextLvl]);
            }

            hierarchy.at(nextLvl).coarseType(X[currLvl], X[nextLvl]);

            clearHistory(nextLvl);
        } else { // postProcess
            if (currLvl == coarsest || hierarchy.getType(currLvl+1) != ELIMINATION) {
                minRes(currLvl, X[currLvl], B[currLvl] - hierarchy.at(currLvl).getLaplacian() * X[currLvl]);
            }

            if (nextLvl > finest) {
                saveIterate(nextLvl, X[nextLvl], B[nextLvl] - hierarchy.at(nextLvl).getLaplacian() * X[nextLvl]);
            }

            if (hierarchy.getType(currLvl) == ELIMINATION) {
                hierarchy.at(currLvl).interpolate(X[currLvl], X[nextLvl], bStages[currLvl]);
            } else {
                Vector xf = X[nextLvl];
                hierarchy.at(currLvl).interpolate(X[currLvl], xf);
                X[nextLvl] += xf;
            }

            if (hierarchy.getType(currLvl) != ELIMINATION) {
                X[nextLvl] = smoother.relax(hierarchy.at(nextLvl).getLaplacian(), B[nextLvl], X[nextLvl], status.numPostSmoothIters);
            }

        }

        currLvl = nextLvl;
    } // while

    // post-cycle finest
    if ((int64_t) hierarchy.size() > finest + 1 && hierarchy.getType(finest+1) != ELIMINATION) { // do an iterate recombination on calculated solutions
        minRes(finest, X[finest], B[finest] - hierarchy.at(finest).getLaplacian() * X[finest]);
    }


    X[finest] -= X[finest].mean();
    x = X[finest];
}

template<class Matrix>
void SolverLamg<Matrix>::saveIterate(index level, const Vector& x, const Vector& r) {
    // update latest pointer
    index i = latestIterate[level];
    latestIterate[level] = (i+1) % MAX_COMBINED_ITERATES;


    // update numIterates
    if (numActiveIterates[level] < MAX_COMBINED_ITERATES) {
        numActiveIterates[level]++;
    }

    // update history array
    history[level][i] = x;
    rHistory[level][i] = r;
}

template<class Matrix>
void SolverLamg<Matrix>::clearHistory(index level) {
    latestIterate[level] = 0;
    numActiveIterates[level] = 0;
}

template<class Matrix>
void SolverLamg<Matrix>::minRes(index level, Vector& x, const Vector& r) {
    if (numActiveIterates[level] > 0) {
        count n = numActiveIterates[level];

        std::vector<index> ARowIdx(r.getDimension()+1);
        std::vector<index> ERowIdx(r.getDimension()+1);

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(r.getDimension()); ++i) {
            for (index k = 0; k < n; ++k) {
                double AEvalue = r[i] - rHistory[level][k][i];
                if (std::fabs(AEvalue) > 1e-25) {
                    ++ARowIdx[i+1];
                }

                double Eval = history[level][k][i] - x[i];
                if (std::fabs(Eval) > 1e-25) {
                    ++ERowIdx[i+1];
                }
            }
        }

        for (index i = 0; i < r.getDimension(); ++i) {
            ARowIdx[i+1] += ARowIdx[i];
            ERowIdx[i+1] += ERowIdx[i];
        }

        std::vector<index> AColumnIdx(ARowIdx[r.getDimension()]);
        std::vector<double> ANonZeros(ARowIdx[r.getDimension()]);

        std::vector<index> EColumnIdx(ERowIdx[r.getDimension()]);
        std::vector<double> ENonZeros(ERowIdx[r.getDimension()]);

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(r.getDimension()); ++i) {
            for (index k = 0, aIdx = ARowIdx[i], eIdx = ERowIdx[i]; k < n; ++k) {
                double AEvalue = r[i] - rHistory[level][k][i];
                if (std::fabs(AEvalue) > 1e-25) {
                    AColumnIdx[aIdx] = k;
                    ANonZeros[aIdx] = AEvalue;
                    ++aIdx;
                }

                double Eval = history[level][k][i] - x[i];
                if (std::fabs(Eval) > 1e-25) {
                    EColumnIdx[eIdx] = k;
                    ENonZeros[eIdx] = Eval;
                    ++eIdx;
                }
            }
        }

        CSRMatrix AE(r.getDimension(), n, ARowIdx, AColumnIdx, ANonZeros, 0.0, true);
        CSRMatrix E(r.getDimension(), n, ERowIdx, EColumnIdx, ENonZeros, 0.0, true);

        Vector alpha = smoother.relax(CSRMatrix::mTmMultiply(AE, AE), CSRMatrix::mTvMultiply(AE, r), Vector(n, 0.0), 10);
        x += E * alpha;
    }

}

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_LAMG_SOLVER_LAMG_HPP_
