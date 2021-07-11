// no-networkit-format
/*
* NeighborhoodFunction.hpp
*
*  Created on: 30.03.2016
*      Author: Maximilian Vogel
*/

#ifndef NETWORKIT_DISTANCE_NEIGHBORHOOD_FUNCTION_HPP_
#define NETWORKIT_DISTANCE_NEIGHBORHOOD_FUNCTION_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 */
class NeighborhoodFunction final : public Algorithm {

public:
    /**
    * Computes the neighborhood function exactly.
    * The neighborhood function N of a graph G for a given distance t is defined
    * as the number of node pairs (u,v) that can be reached within distance t.
    *
    * @param G the given graph
    * @return the exact effective diameter of the graph
    */
    NeighborhoodFunction(const Graph& G);

    void run() override;

    /**
     * Returns the neighborhood function of the graph.
     * @return the neighborhood function of the graph
     */
    std::vector<count> getNeighborhoodFunction() const;

private:
    const Graph* G;
    std::vector<count> result;
};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_NEIGHBORHOOD_FUNCTION_HPP_
