// no-networkit-format
/*
 * Eccentricity.hpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#ifndef NETWORKIT_DISTANCE_ECCENTRICITY_HPP_
#define NETWORKIT_DISTANCE_ECCENTRICITY_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * The eccentricity of a node `u` is defined as the distance to the farthest
 * node from `u`. In other words, it is the longest shortest-path starting from
 * node `u`.
 */
class Eccentricity {

public:

    /**
     * @return The farthest node v, and the length of the shortest path to v.
     */
    static std::pair<node, count> getValue(const Graph& G, node u);
};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_ECCENTRICITY_HPP_
