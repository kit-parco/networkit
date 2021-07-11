/*
 * AlgebraicGlobals.h
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_ALGEBRAIC_ALGEBRAIC_GLOBALS_HPP_
#define NETWORKIT_ALGEBRAIC_ALGEBRAIC_GLOBALS_HPP_

#include <networkit/Globals.hpp>

namespace NetworKit {

/** Represents a matrix entry s.t. matrix(row, column) = value */
struct Triplet {
    index row;
    index column;
    double value;
};

/** Floating point epsilon to use in comparisons. */
constexpr double FLOAT_EPSILON = 1e-9;

} // namespace NetworKit

#endif // NETWORKIT_ALGEBRAIC_ALGEBRAIC_GLOBALS_HPP_
