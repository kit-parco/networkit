// no-networkit-format
/*
 * LAMGSettings.h
 *
 *  Created on: 09.01.2015
 *      Author: Michael
 */

#ifndef NETWORKIT_NUMERICS_LAMG_LAMG_SETTINGS_HPP_
#define NETWORKIT_NUMERICS_LAMG_LAMG_SETTINGS_HPP_

#include <networkit/Globals.hpp>

namespace NetworKit {

// #initial test vectors (TVs) at each level.
constexpr count TV_NUM = 4;
// #TVs to add upon each aggregation coarsening
constexpr count TV_INC = 1;
// maximum allowed#TVs
constexpr count TV_MAX = 10;
// #global sweeps to perform on each initial TV
constexpr count SETUP_TV_SWEEPS = 4;
// max size for direct solver
constexpr count MAX_DIRECT_SOLVE_SIZE = 200;
// maximum #aggregation coarsening levels to construct during the setup phase
constexpr count SETUP_MAX_AGG_LEVELS = 100;
// maximum TOTAL # of coarsening levels to construct during the setup phase
constexpr count SETUP_MAX_LEVELS = 100;
// solution cycle index. Also the design cycle index during setup phase.
constexpr double SETUP_CYCLE_INDEX = 1.5;
// minimum number of sweeps to run to estimate relax ACF
constexpr count SETUP_RELAX_ACF_MIN_SWEEPS = 7;
// if relaxation converges at this rate or faster, this becomes the CoarsestLevel
constexpr double SETUP_MAX_COARSE_RELAX_ACF = 0.3;

constexpr count MAX_COMBINED_ITERATES = 4;

/**************************
 * SETUP - Elimination    *
 **************************/

// maximum degree up to which a node gets eliminated
constexpr count SETUP_ELIMINATION_MAX_DEGREE = 4;
// maximum number of elimination stages
constexpr count SETUP_ELIMINATION_MAX_STAGES = 1000;
// node elimination stops if number of coarse nodes is only this fraction of the initial number of nodes
constexpr double SETUP_ELIMINATION_MIN_ELIM_FRACTION = 0.01;

/**************************
 * SETUP - Aggregation    *
 **************************/

constexpr double SETUP_AGGREGATION_WEAK_EDGE_THRESHOLD = 0.1;
// all nodes with a degree equal or higher than this threshold are marked as high-degree nodes
constexpr count SETUP_AGGREGATION_DEGREE_THRESHOLD = 8;
// #sweeps (nu) to use during coarsening stages
constexpr count SETUP_NU_DEFAULT = 3;

constexpr double SETUP_COARSENING_WORK_GUARD = 0.7;

constexpr count SETUP_MIN_AGGREGATION_STAGES = 1;

constexpr count SETUP_MAX_AGGREGATION_STAGES = 2;

constexpr count SETUP_RELAX_COARSEST_SWEEPS = 400;

}



#endif // NETWORKIT_NUMERICS_LAMG_LAMG_SETTINGS_HPP_
