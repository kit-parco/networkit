// no-networkit-format
/*
 * LevelAggregation.h
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#ifndef NETWORKIT_NUMERICS_LAMG_LEVEL_LEVEL_AGGREGATION_HPP_
#define NETWORKIT_NUMERICS_LAMG_LEVEL_LEVEL_AGGREGATION_HPP_

#include <networkit/numerics/LAMG/Level/Level.hpp>

namespace NetworKit {

/**
 * @ingroup numerics
 */
template<class Matrix>
class LevelAggregation : public Level<Matrix> {
private:
    Matrix P; // interpolation matrix (n x nc)
    Matrix R; // restriction matrix (nc x n)

public:
    LevelAggregation(const Matrix& A, const Matrix& P, const Matrix& R) : Level<Matrix>(LevelType::AGGREGATION, A), P(P), R(R) {}

    void coarseType(const Vector &xf, Vector &xc) const override;

    void restrict(const Vector &bf, Vector &bc) const override;

    void interpolate(const Vector &xc, Vector &xf) const override;
};

template<class Matrix>
void LevelAggregation<Matrix>::coarseType(const Vector& /*xf*/, Vector& xc) const {
    xc = Vector(P.numberOfColumns(), 0.0);
}

template<class Matrix>
void LevelAggregation<Matrix>::restrict(const Vector& bf, Vector& bc) const {
    bc = R * bf;
}

template<class Matrix>
void LevelAggregation<Matrix>::interpolate(const Vector& xc, Vector& xf) const {
    xf = P * xc;
}

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_LAMG_LEVEL_LEVEL_AGGREGATION_HPP_
