// no-networkit-format
/*
 * EpidemicSimulationSEIR.cpp
 *
 *  Created on: 20.11.2015
 *      Author: Christian Staudt
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/simulation/EpidemicSimulationSEIR.hpp>

namespace NetworKit {

EpidemicSimulationSEIR::EpidemicSimulationSEIR(const Graph& G, count tMax, double transP, count eTime, count iTime, node zero) : Algorithm(), G(&G), tMax(tMax), transP(transP), eTime(eTime), iTime(iTime), zero(zero)  {}

void EpidemicSimulationSEIR::run() {

    using State = EpidemicSimulationSEIR::State;

    index t = 0;

    //initialize state and timestamp arrays
    state.resize(G->upperNodeIdBound(), State::U);
    timestamp.resize(G->upperNodeIdBound(), none);

    auto setState = [&](node v, State X){
        state[v] = X;
        timestamp[v] = t;
    };

    // initialize nodes to Susceptible
    G->parallelForNodes([&](node v) {
        setState(v, State::S);
    });

    // contact may expose susceptible node to infection
    auto contact = [&](node v) {
        if ((state[v] == State::S) && (Aux::Random::probability() <= transP)) {
            setState(v, State::E);
        }
    };

    // update state of nodes
    auto sweep = [&](node u) {
        if (state[u] == State::S) {
            // do nothing
        } else if (state[u] == State::E) {
            // exposed nodes become infectious after time
            if ((t - timestamp[u]) >= eTime) {
                setState(u, State::I);
            }
        } else if (state[u] == State::I) {
            // contact neighbors of infectious node
            G->forNeighborsOf(u, [&](node v){
                contact(v);
            });
            // infectious nodes become removed after time
            if ((t - timestamp[u]) >= iTime) {
                setState(u, State::R);
            }
        } else if (state[u] == State::R) {
            // do nothing
        } else if (state[u] == State::U) {
            throw std::runtime_error("node in undefined state encountered - should not happen");
        } else {
            throw std::runtime_error("else branch taken - should not happen");
        }
    };


    auto census = [&]() {
        std::vector<count> data(5);
        G->forNodes([&](node v) {
            data[(index) state[v]] += 1;
        });
        return data;
    };


    // if starting node node provided, start with random node
    if (zero == none) {
        zero = GraphTools::randomNode(*G);
    }
    INFO("zero node: ", zero);
    setState(zero, State::I);	// infect node zero

    while (t < tMax) {
        G->parallelForNodes(sweep);
        auto populations = census();

        for (int s = (int) State::S; s != (int) State::U; ++s) {
            std::vector<count> data = {zero, t, (count)s, populations[s]};
            stats.push_back(data);
        }

        t += 1;
    }

    hasRun = true;
}


std::vector<std::vector<count>> EpidemicSimulationSEIR::getData() const {
    return stats;
}

} /* namespace NetworKit */
