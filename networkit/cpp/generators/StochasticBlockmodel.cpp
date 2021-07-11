// no-networkit-format
/*
* StochasticBlockmodel.cpp
*
*  Created on: 13.08.2014
*      Author: Christian Staudt
*/

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/StochasticBlockmodel.hpp>

namespace NetworKit {

StochasticBlockmodel::StochasticBlockmodel(count n, count nBlocks, const std::vector<index>& membership, const std::vector<std::vector<double>>& affinity)
    : n(n), membership(membership), affinity(affinity) {

    std::string errorMessage = "affinity matrix must be of size nBlocks x nBlocks";
    if (affinity.size() != nBlocks) {
        throw std::runtime_error(errorMessage);
    }

    for (const auto &row : affinity) {
        if (row.size() != nBlocks) {
            throw std::runtime_error(errorMessage);
        }
    }

    if (membership.size() != n) {
        throw std::runtime_error("membership list must be of size nNodes");
    }
}


Graph StochasticBlockmodel::generate() {
    Graph G(n);

    G.forNodePairs([&](node u, node v) {
        index a = membership.at(u);
        index b = membership.at(v);
        assert (a < affinity.size());
        assert (b < affinity.size());
        double p = affinity.at(a).at(b);
        double r = Aux::Random::real();
        if (r <= p) {
            G.addEdge(u, v);
        }
    });
    return G;
}

} /* namespace NetworKit */
