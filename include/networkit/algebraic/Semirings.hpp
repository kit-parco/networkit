// no-networkit-format
/*
 * Semirings.hpp
 *
 *  Created on: May 31, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_ALGEBRAIC_SEMIRINGS_HPP_
#define NETWORKIT_ALGEBRAIC_SEMIRINGS_HPP_

#include <algorithm>

// *****************************************************
// 					Semiring Definitions
// *****************************************************

/**
 * @ingroup algebraic
 * add: arithmetic add
 * mult: arithmetic multiplication
 * zero: 0
 * one: 1
 * codomain = (-infty, +infty)
 */
class ArithmeticSemiring final {
public:
    ArithmeticSemiring() = default;
    ~ArithmeticSemiring() = default;

    inline static double add(double a, double b) {
        return a + b;
    }

    inline static double mult(double a, double b) {
        return a * b;
    }

    inline static double zero() {return 0;};

    inline static double one() {return 1;};
};

/**
 * @ingroup algebraic
 * add: min
 * mult: arithmetic add
 * zero: +infty
 * one: 0
 * codomain = (-infty, +infty]
 */
class MinPlusSemiring final {
public:
    MinPlusSemiring() = default;
    ~MinPlusSemiring() = default;

    inline static double add(double a, double b) {
        return std::min(a,b);
    }

    inline static double mult(double a, double b) {
        return a+b;
    }

    inline static double zero() {return std::numeric_limits<double>::infinity();};

    inline static double one() {return 0;};
};

/**
 * @ingroup algebraic
 * add: max
 * mult: arithmetic add
 * zero: -infty
 * one: 0
 * codomain = [-infty, +infty)
 */
class MaxPlusSemiring final {
public:
    MaxPlusSemiring() = default;
    ~MaxPlusSemiring() = default;

    inline static double add(double a, double b) {
        return std::max(a,b);
    }

    inline static double mult(double a, double b) {
        return a+b;
    }

    inline static double zero() {return -std::numeric_limits<double>::infinity();};

    inline static double one() {return 0;};
};

/**
 * @ingroup algebraic
 * add: min
 * mult: max
 * zero: +infty
 * one: -infty
 * codomain = [-infty, +infty]
 */
class MinMaxSemiring final {
public:
    MinMaxSemiring() = default;
    ~MinMaxSemiring() = default;

    inline static double add(double a, double b) {
        return std::min(a,b);
    }

    inline static double mult(double a, double b) {
        return std::max(a,b);
    }

    inline static double zero() {return std::numeric_limits<double>::infinity();};

    inline static double one() {return -std::numeric_limits<double>::infinity();};
};

/**
 * @ingroup algebraic
 * add: max
 * mult: min
 * zero: -infty
 * one: +infty
 * codomain = [-infty, +infty]
 */
class MaxMinSemiring final {
public:
    MaxMinSemiring() = default;
    ~MaxMinSemiring() = default;

    inline static double add(double a, double b) {
        return std::max(a,b);
    }

    inline static double mult(double a, double b) {
        return std::min(a,b);
    }

    inline static double zero() {return -std::numeric_limits<double>::infinity();};

    inline static double one() {return 0;};
};

/**
 * @ingroup algebraic
 * add: logical or
 * mult: logical and
 * zero: 0
 * one: 1
 * codomain = [-infty, +infty]
 */
class IntLogicalSemiring final {
public:
    IntLogicalSemiring() = default;
    ~IntLogicalSemiring() = default;

    inline static double add(double a, double b) {
        return (int) a || (int) b;
    }

    inline static double mult(double a, double b) {
        return (int) a && (int) b;
    }

    inline static double zero() {return 0;};

    inline static double one() {return 1;};
};

/**
 * @ingroup algebraic
 * add: xor
 * mult: bitwise and
 * zero: 0
 * one: 1
 * codomain = [0, 1]
 */
class GaloisFieldSemiring final {
public:
    GaloisFieldSemiring() = default;
    ~GaloisFieldSemiring() = default;

    inline static double add(double a, double b) {
        return (int) a ^ (int) b;
    }

    inline static double mult(double a, double b) {
        return (int) a & (int) b;
    }

    inline static double zero() {return 0;};

    inline static double one() {return 1;};
};



#endif // NETWORKIT_ALGEBRAIC_SEMIRINGS_HPP_
