// no-networkit-format
/*
 * RandomNodeEdgeScore.hpp
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#ifndef NETWORKIT_SPARSIFICATION_RANDOM_NODE_EDGE_SCORE_HPP_
#define NETWORKIT_SPARSIFICATION_RANDOM_NODE_EDGE_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

class RandomNodeEdgeScore final : public EdgeScore<double> {

public:
    RandomNodeEdgeScore(const Graph& graph, double rneRatio = 0.8);
    void run() override;
    double score(edgeid eid) override;
    double score(node u, node v) override;

private:
    double rneRatio;

};

} /* namespace NetworKit */

#endif // NETWORKIT_SPARSIFICATION_RANDOM_NODE_EDGE_SCORE_HPP_
