// no-networkit-format
/*
*  HopPlotApproximation.hpp
*
*  Created on: 30.03.2016
*      Author: Maximilian Vogel
*/

#ifndef NETWORKIT_DISTANCE_HOP_PLOT_APPROXIMATION_HPP_
#define NETWORKIT_DISTANCE_HOP_PLOT_APPROXIMATION_HPP_

#include <map>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 */
class HopPlotApproximation final : public Algorithm {

public:
    /**
    * Computes an approximation of the hop-plot of a given graph
    * The hop-plot is the set of pairs (d, g(g)) for each natural number d
    * and where g(d) is the fraction of connected node pairs whose shortest connecting path has length at most d.
    * Implementation after the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining in Massive Graphs"[1]
    *
    * [1] by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf
    *
    * @param G the given graph
    * @param maxDistance the maximum path length that shall be considered. set 0 for infinity/diameter of the graph
    * @param k the number of parallel approximations to get a more robust result; default = 64
    * @param r the amount of bits that are added to the length of the bitmask to improve the accuracy; default = 7
    * @return the approximated hop-plot of the graph
    */
    HopPlotApproximation(const Graph& G, count maxDistance = 0, count k = 64, count r = 7);

    void run() override;

    /**
     * Returns the approximated hop-plot of the graph.
     * @return the approximated hop-plot of the graph
     */
    std::map<count, double> getHopPlot() const;

private:
    const Graph* G;
    const count maxDistance;
    const count k;
    const count r;
    std::map<count, double> hopPlot;

};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_HOP_PLOT_APPROXIMATION_HPP_
