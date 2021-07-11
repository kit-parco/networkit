// no-networkit-format
/*
 * MissingMath.h
 *
 *  Created on: 21.03.2013
 *      Author: cls
 */

#ifndef NETWORKIT_AUXILIARY_MISSING_MATH_HPP_
#define NETWORKIT_AUXILIARY_MISSING_MATH_HPP_

#include <cinttypes>
#include <cmath>
#include <cassert>
#include <stdexcept>

namespace Aux {

/**
 * Math functions not provided by the standard library.
 */
namespace MissingMath {


inline int64_t binomial(int64_t n, int64_t k) {
    if (k == 0) return 1;
    if (2 * k > n) return binomial(n, n - k);

    int64_t b = n - k + 1;
    for (int64_t i = 2; i <= k; ++i) {
        b = b * (n - k + i);
        b = b / i;
    }
    return b;
}

inline double log_b(double x, double b) {
    if (x == 0) {
        throw std::domain_error("log(0) is undefined");
    }
    assert (log(b) != 0);
    return log(x) / log(b);
}

} /* namespace MissingMath */

} /* namespace Aux */
#endif // NETWORKIT_AUXILIARY_MISSING_MATH_HPP_
