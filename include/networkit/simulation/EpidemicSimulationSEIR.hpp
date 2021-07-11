// no-networkit-format
/*
 * EpidemicSimulationSEIR.hpp
 *
 *  Created on: 20.11.2015
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_SIMULATION_EPIDEMIC_SIMULATION_SEIR_HPP_
#define NETWORKIT_SIMULATION_EPIDEMIC_SIMULATION_SEIR_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup simulation
 * Node centrality index which ranks nodes by their degree.
 * Optional normalization by maximum degree.
 */
class EpidemicSimulationSEIR final : public Algorithm {
public:
    /**
     * @param G The network.
     */
    EpidemicSimulationSEIR(const Graph& G, count tMax, double transP, count eTime, count iTime, node zero);

    void run() override;

    std::vector<std::vector<count>> getData() const;

private:

    const Graph* G;
    count tMax;
    double transP;
    count eTime;
    count iTime;
    node zero;
    enum class State {S, E, I, R, U}; // Susceptible, Exposed, Infectious, Removed, Undefined
    std::vector<State> state;
    std::vector<index> timestamp;
    std::vector<std::vector<count>> stats;
};

} /* namespace NetworKit */

#endif // NETWORKIT_SIMULATION_EPIDEMIC_SIMULATION_SEIR_HPP_
