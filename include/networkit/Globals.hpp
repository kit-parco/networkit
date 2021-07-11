// no-networkit-format
/*
 * Globals.h
 *
 *  Created on: 06.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_GLOBALS_HPP_
#define NETWORKIT_GLOBALS_HPP_

#include <cstdint>
#include <limits>
#include <utility>

namespace NetworKit {
    using index = uint64_t; ///< more expressive name for an index into an array

    /// Should be used in OpenMP parallel for-loops and is associated with unsigned semantics.
    /// On MSVC it falls back to being signed, as MSVC does not support unsigned parallel fors.
#ifdef _MSC_VER
    using omp_index = int64_t;
#else
    using omp_index = index;
#endif // _MSC_VER

    using count      = uint64_t; ///< more expressive name for an integer quantity
    using node       = index;    ///< node indices are 0-based
    using edgeweight = double;   ///< edge weight type
    using edgeid     = index;    ///< edge id

    constexpr index none = std::numeric_limits<index>::max(); ///< value for not existing nodes/edges
    constexpr edgeweight defaultEdgeWeight = 1.0;
    constexpr edgeweight nullWeight = 0.0;

    using coordinate = double;

    constexpr double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

}

#endif // NETWORKIT_GLOBALS_HPP_
