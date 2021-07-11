// no-networkit-format
/*
* NeighborhoodFunctionHeuristic.hpp
*
*      Author: Maximilian Vogel
*/

#ifndef NETWORKIT_DISTANCE_NEIGHBORHOOD_FUNCTION_HEURISTIC_HPP_
#define NETWORKIT_DISTANCE_NEIGHBORHOOD_FUNCTION_HEURISTIC_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 */
class NeighborhoodFunctionHeuristic final : public Algorithm {

public:
    enum SelectionStrategy {
        RANDOM,
        SPLIT
    };
    /**
    * Computes a heuristic of the neighborhood function.
    *
    * The algorithm runs nSamples breadth-first searches and scales the results up to the actual amount of nodes.
    * Accepted strategies are "split" and "random".
    *
    * @param G the given graph
    * @param nSamples the amount of samples, set to zero for heuristic of max(sqrt(m), 0.15*n)
    * @param strategy the strategy to select the samples, accepts "random" or "split"
    */
    NeighborhoodFunctionHeuristic(const Graph& G, count nSamples = 0, SelectionStrategy strategy = SPLIT);

    void run() override;

    /**
     * Returns the approximated neighborhood function of the graph.
     * @return the approximated neighborhood function of the graph
     */
    std::vector<count> getNeighborhoodFunction() const;

private:
    const Graph* G;
    const count nSamples;
    const SelectionStrategy strategy;
    std::vector<count> result;

    /* selection schemes implemented as private functions */
    std::vector<node> random(const Graph& G, count nSamples);
    std::vector<node> split(const Graph& G, count nSamples);

};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_NEIGHBORHOOD_FUNCTION_HEURISTIC_HPP_
