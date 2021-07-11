// no-networkit-format
/*
 * Smoother.h
 *
 *  Created on: 31.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_NUMERICS_SMOOTHER_HPP_
#define NETWORKIT_NUMERICS_SMOOTHER_HPP_

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>

#include <limits>
#include <networkit/algebraic/DynamicMatrix.hpp>

namespace NetworKit {

/**
 * @ingroup numerics
 * Abstract base class of a smoother.
 */
template<class Matrix>
class Smoother {
public:
    Smoother() {}
    virtual ~Smoother(){}

    virtual Vector relax(const Matrix &A, const Vector &b, const Vector &initialGuess,
                         count maxIterations = std::numeric_limits<count>::max()) const = 0;
    virtual Vector relax(const Matrix &A, const Vector &b,
                         count maxIterations = std::numeric_limits<count>::max()) const = 0;
};

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_SMOOTHER_HPP_
