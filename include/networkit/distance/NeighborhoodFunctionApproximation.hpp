// no-networkit-format
/*
* NeighborhoodFunctionApproximation.hpp
*
*  Created on: 30.03.2016
*      Author: Maximilian Vogel
*/

#ifndef NETWORKIT_DISTANCE_NEIGHBORHOOD_FUNCTION_APPROXIMATION_HPP_
#define NETWORKIT_DISTANCE_NEIGHBORHOOD_FUNCTION_APPROXIMATION_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 */
class NeighborhoodFunctionApproximation final : public Algorithm {

public:
    /**
    * Computes an approximation of the neighborhood function.
    * The neighborhood function N of a graph G for a given distance t is defined
    * as the number of node pairs (u,v) that can be reached within distance t.
    * Implementation after the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining in Massive Graphs"[1]
    *
    * [1] by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf
    *
    * @param G the given graph
    * @param k the number of parallel approximations to get a more robust result; default = 64
    * @param r the amount of bits that are added to the length of the bitmask to improve the accuracy; default = 7
    */
    NeighborhoodFunctionApproximation(const Graph& G, count k = 64, count r = 7);

    void run() override;

    /**
     * Returns the approximated neighborhood function of the graph.
     * @return the approximated neighborhood function of the graph
     */
    std::vector<count> getNeighborhoodFunction() const;

private:
    const Graph* G;
    const count k;
    const count r;
    std::vector<count> result;
};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_NEIGHBORHOOD_FUNCTION_APPROXIMATION_HPP_
