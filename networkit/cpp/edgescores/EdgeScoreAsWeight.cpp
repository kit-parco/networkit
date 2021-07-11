// no-networkit-format
/*
 * EdgeScoreAsWeight.cpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#include <networkit/edgescores/EdgeScoreAsWeight.hpp>

namespace NetworKit {

EdgeScoreAsWeight::EdgeScoreAsWeight(const Graph& G, const std::vector<double>& score, bool squared, edgeweight offset, edgeweight factor) :
        G(&G), score(&score), squared(squared), offset(offset), factor(factor) {
}

Graph EdgeScoreAsWeight::calculate() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    Graph result(*G, true, false);

    if (squared) {
        G->parallelForEdges([&](node u, node v, edgeid eid) {
            result.setWeight(u, v, offset + factor * (*score)[eid] * (*score)[eid]);
        });
    } else {
        G->parallelForEdges([&](node u, node v, edgeid eid) {
            result.setWeight(u, v, offset + factor * (*score)[eid]);
        });
    }

    return result;
}

} // namespace NetworKit
