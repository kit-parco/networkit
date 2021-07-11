// no-networkit-format
/*
*  EffectiveDiameter.hpp
*
*  Created on: 16.06.2014
*      Author: Marc Nemes
*/

#ifndef NETWORKIT_DISTANCE_EFFECTIVE_DIAMETER_HPP_
#define NETWORKIT_DISTANCE_EFFECTIVE_DIAMETER_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 */
class EffectiveDiameter final : public Algorithm {

public:
    /**
    * Computes the effective diameter exactly.
    * The effective diameter is defined as the number of edges on average to reach \p ratio of all other nodes.
    * @param G the given graph
    * @param ratio the ratio of nodes that should be connected (0,1], default = 0.9
    */
    EffectiveDiameter(const Graph& G, double ratio = 0.9);

    void run() override;

    /**
     * Returns the exact effective diameter of the graph.
     * @return the exact effective diameter of the graph
     */
    double getEffectiveDiameter() const;

private:
    const Graph* G;
    const double ratio;
    double effectiveDiameter;
};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_EFFECTIVE_DIAMETER_HPP_
