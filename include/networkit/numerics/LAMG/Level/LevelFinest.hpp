// no-networkit-format
/*
 * LevelFinest.h
 *
 *  Created on: 10.01.2015
 *      Author: Michael
 */

#ifndef NETWORKIT_NUMERICS_LAMG_LEVEL_LEVEL_FINEST_HPP_
#define NETWORKIT_NUMERICS_LAMG_LEVEL_LEVEL_FINEST_HPP_

#include <networkit/numerics/LAMG/Level/Level.hpp>

namespace NetworKit {

/**
 * @ingroup numerics
 */
template<class Matrix>
class LevelFinest : public Level<Matrix> {
public:
    LevelFinest() : Level<Matrix>(LevelType::FINEST) {}
    LevelFinest(const Matrix& A) : Level<Matrix>(LevelType::FINEST, A) {}
    ~LevelFinest() override = default;
};

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_LAMG_LEVEL_LEVEL_FINEST_HPP_
