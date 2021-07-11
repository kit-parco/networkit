// no-networkit-format
/*
 * EliminationStage.h
 *
 *  Created on: 09.01.2015
 *      Author: Michael
 */

#ifndef NETWORKIT_NUMERICS_LAMG_LEVEL_ELIMINATION_STAGE_HPP_
#define NETWORKIT_NUMERICS_LAMG_LEVEL_ELIMINATION_STAGE_HPP_

#include <networkit/algebraic/Vector.hpp>

namespace NetworKit {

/**
 * @ingroup numerics
 */
template<class Matrix>
class EliminationStage {
private:
    Matrix P; // interpolation matrix
    Matrix R;
    Vector q; // coarse result correction vector
    std::vector<index> fSet;
    std::vector<index> cSet;

public:
    EliminationStage(const Matrix& P, const Vector& q, const std::vector<index>& fSet, const std::vector<index>& cSet) : P(P), R(P.transpose()), q(q), fSet(fSet), cSet(cSet) {}

    inline const Matrix& getP() const {
        return P;
    }

    inline const Matrix& getR() const {
        return R;
    }

    inline const Vector& getQ() const {
        return q;
    }

    inline const std::vector<index>& getFSet() const {
        return fSet;
    }

    inline const std::vector<index>& getCSet() const {
        return cSet;
    }

    inline count getN() const {
        return fSet.size() + cSet.size();
    }
};

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_LAMG_LEVEL_ELIMINATION_STAGE_HPP_
