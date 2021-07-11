// no-networkit-format
/*
 * EdgeScoreAsWeight.hpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef NETWORKIT_EDGESCORES_EDGE_SCORE_AS_WEIGHT_HPP_
#define NETWORKIT_EDGESCORES_EDGE_SCORE_AS_WEIGHT_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class EdgeScoreAsWeight final {

public:
    EdgeScoreAsWeight(const Graph& G, const std::vector<double>& score, bool squared = false, edgeweight offset = 1, edgeweight factor = 1);
    Graph calculate();

private:
    const Graph* G;
    const std::vector<double>* score;
    bool squared;
    edgeweight offset;
    edgeweight factor;

};

} // namespace NetworKit

#endif // NETWORKIT_EDGESCORES_EDGE_SCORE_AS_WEIGHT_HPP_
