// no-networkit-format
/*
 * GraphStructuralRandMeasure.cpp
 *
 *  Created on: 16.04.2013
 *      Author: cls
 */

#include <networkit/community/GraphStructuralRandMeasure.hpp>

namespace NetworKit {

double GraphStructuralRandMeasure::getDissimilarity(const Graph& G, const Partition& first,
        const Partition& second) {
    count m = G.numberOfEdges();
    if (m == 0) {
        throw std::runtime_error("The graph-structural rand measure is not defined for graphs without edges.");
    }

    double e11 = 0.0; // number of connected node pairs for which clusterings agree
    double e00 = 0.0; // number of connected node pairs for which clusterings disagree

    G.forEdges([&](node u, node v){
        if ((first[u] == first[v]) && (second[u] == second[v])) {
            e11 += 1;
        } else if ((first[u] != first[v]) && (second[u] != second[v])) {
            e00 += 1;
        }
    });

    double rand = 1.0 - (e11 + e00) * (1.0 / m);

    // assert range [0, 1]
    assert (rand <= 1.0);
    assert (rand >= 0.0);
    return rand;
}

} /* namespace NetworKit */
