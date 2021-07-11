// no-networkit-format
/*
 * LevelElimination.h
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#ifndef NETWORKIT_NUMERICS_LAMG_LEVEL_LEVEL_ELIMINATION_HPP_
#define NETWORKIT_NUMERICS_LAMG_LEVEL_LEVEL_ELIMINATION_HPP_

#include <networkit/numerics/LAMG/Level/Level.hpp>
#include <networkit/numerics/LAMG/Level/EliminationStage.hpp>

namespace NetworKit {

/**
 * @ingroup numerics
 */
template<class Matrix>
class LevelElimination : public Level<Matrix> {
private:
    std::vector<EliminationStage<Matrix>> coarseningStages;
    std::vector<index> cIndexFine;

    void subVectorExtract(Vector& subVector, const Vector& vector, const std::vector<index>& elements) const;

public:
    LevelElimination(const Matrix& A, const std::vector<EliminationStage<Matrix>>& coarseningStages);

    void coarseType(const Vector &xf, Vector &xc) const override;
    void restrict(const Vector &bf, Vector &bc, std::vector<Vector> &bStages) const override;
    void interpolate(const Vector &xc, Vector &xf,
                     const std::vector<Vector> &bStages) const override;
};

template<class Matrix>
LevelElimination<Matrix>::LevelElimination(const Matrix& A, const std::vector<EliminationStage<Matrix>>& coarseningStages) : Level<Matrix>(LevelType::ELIMINATION, A), coarseningStages(coarseningStages) {
    cIndexFine = std::vector<index>(this->A.numberOfRows());
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(cIndexFine.size()); ++i) {
        cIndexFine[i] = i;
    }

    for (index k = coarseningStages.size(); k-- > 0;) {
        for (index i = 0; i < cIndexFine.size(); ++i) {
            assert(cIndexFine[i] < coarseningStages[k].getCSet().size());
            cIndexFine[i] = coarseningStages[k].getCSet()[cIndexFine[i]];
        }
    }
}

template<class Matrix>
void LevelElimination<Matrix>::coarseType(const Vector& xf, Vector& xc) const {
    xc = Vector(this->A.numberOfRows());
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(xc.getDimension()); ++i) {
        xc[i] = xf[cIndexFine[i]];
    }
}

template<class Matrix>
void LevelElimination<Matrix>::restrict(const Vector& bf, Vector& bc, std::vector<Vector>& bStages) const {
    bStages.resize(coarseningStages.size() + 1);
    bStages[0] = bf;
    bc = bf;
    index curStage = 0;
    for (const EliminationStage<Matrix>& s : coarseningStages) {
        //Vector bOld = bStages[curStage];
        Vector bCSet;
        subVectorExtract(bCSet, bc, s.getCSet());

        Vector bFSet;
        subVectorExtract(bFSet, bc, s.getFSet());
        bc = bCSet + s.getR() * bFSet;
        bStages[curStage+1] = bc; // b = b.c + s.P^T * b.f

        curStage++;
    }
}

template<class Matrix>
void LevelElimination<Matrix>::interpolate(const Vector& xc, Vector& xf, const std::vector<Vector>& bStages) const {
    Vector currX = xc;
    for (index k = coarseningStages.size(); k-- > 0;) {
        const EliminationStage<Matrix>& s = coarseningStages[k];
        xf = Vector(s.getN());
        Vector bFSet;
        subVectorExtract(bFSet, bStages[k], s.getFSet());

        Vector bq(bFSet.getDimension());
        const Vector &q = s.getQ();
#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(bq.getDimension()); ++i) { // bq = s.q .* b.f
            bq[i] = q[i] * bFSet[i];
        }
        Vector xFSet = s.getP() * currX + bq;

        const std::vector<index> &fSet = s.getFSet();
#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(xFSet.getDimension()); ++i) {
            xf[fSet[i]] = xFSet[i];
        }

        const std::vector<index> &cSet = s.getCSet();
#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(currX.getDimension()); ++i) {
            xf[cSet[i]] = currX[i];
        }

        currX = xf;
    }
}

template<class Matrix>
void LevelElimination<Matrix>::subVectorExtract(Vector& subVector, const Vector& vector, const std::vector<index>& elements) const {
    subVector = Vector(elements.size());
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(elements.size()); ++i) {
        subVector[i] = vector[elements[i]];
    }
}

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_LAMG_LEVEL_LEVEL_ELIMINATION_HPP_
