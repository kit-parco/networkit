// no-networkit-format
/*
 * Multiprecision.hpp
 *
 *  Created on: 21.03.2013
 *      Author: cls
 */

#ifndef NETWORKIT_AUXILIARY_MULTIPRECISION_HPP_
#define NETWORKIT_AUXILIARY_MULTIPRECISION_HPP_

#include <ttmath/ttmath.hpp>

namespace NetworKit {
    using bigfloat = ttmath::Big<TTMATH_BITS(64),TTMATH_BITS(64)>;	///< big floating point number
}

#endif // NETWORKIT_AUXILIARY_MULTIPRECISION_HPP_
