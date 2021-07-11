// no-networkit-format
/*
 * Matcher.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt
 */

#include <networkit/matching/Matcher.hpp>

namespace NetworKit {

Matcher::Matcher(const Graph& G): G(&G), M(G.upperNodeIdBound()), edgeScoresAsWeights(false) {}

Matcher::Matcher(const Graph& G, const std::vector<double>& edgeScores): G(&G), M(G.upperNodeIdBound()), edgeScoresAsWeights(true), edgeScores(edgeScores) {
    if (!G.hasEdgeIds()) throw std::invalid_argument("index edges of input graph first");
}

Matching Matcher::getMatching() const {
    return M;
}

} /* namespace NetworKit */
