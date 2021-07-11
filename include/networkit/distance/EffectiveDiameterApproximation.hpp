// no-networkit-format
/*
*  EffectiveDiameterApproximation.hpp
*
*  Created on: 29.03.16
*      Author: Maximilian Vogel
*/

#ifndef NETWORKIT_DISTANCE_EFFECTIVE_DIAMETER_APPROXIMATION_HPP_
#define NETWORKIT_DISTANCE_EFFECTIVE_DIAMETER_APPROXIMATION_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 */
class EffectiveDiameterApproximation final : public Algorithm {

public:
    /**
    * Approximates the effective diameter of a given graph.
    * The effective diameter is defined as the number of edges on average to reach \p ratio of all other nodes.
    * Implementation after the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining in Massive Graphs"[1]
    *
    * [1] by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf
    *
    * @param G the given graph
    * @param ratio the ratio of nodes that should be connected (0,1]; default = 0.9
    * @param k the number of parallel approximations to get a more robust result; default = 64
    * @param r the amount of bits that are added to the length of the bitmask to improve the accuracy; default = 7
    */
    EffectiveDiameterApproximation(const Graph &G, double ratio = 0.9, count k = 64, count r = 7);

    void run() override;

    /**
     * Returns the exact effective diameter of the graph.
     * @return the exact effective diameter of the graph
     */
    double getEffectiveDiameter() const;

private:
    const Graph* G;
    const double ratio;
    const count k;
    const count r;
    double effectiveDiameter;
};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_EFFECTIVE_DIAMETER_APPROXIMATION_HPP_
