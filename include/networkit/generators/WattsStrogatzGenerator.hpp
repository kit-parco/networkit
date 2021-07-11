/*
 * WattsStrogatzGenerator.hpp
 *
 *  Created on: 09.07.2014
 *      Author: Simon Bischof
 */
#ifndef NETWORKIT_GENERATORS_WATTS_STROGATZ_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_WATTS_STROGATZ_GENERATOR_HPP_

#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 */
class WattsStrogatzGenerator final : public StaticGraphGenerator {

public:
    /**
     * Constructs a graph according to the Watts and Strogatz model
     * (https://en.wikipedia.org/wiki/Watts_and_Strogatz_model),
     * which produces graphs with high clustering and low average path length.
     *
     * First, a regular ring lattice is generated.
     * Then some edges are rewired randomly.
     *
     * @param nNodes number of nodes in target graph
     * @param nNeighbors number of neighbors on each side of a node
     * @param p rewiring probability
     */
    WattsStrogatzGenerator(count nNodes, count nNeighbors, double p);

    Graph generate() override;

private:
    count nNodes;
    count nNeighbors;
    double p;
};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_WATTS_STROGATZ_GENERATOR_HPP_
