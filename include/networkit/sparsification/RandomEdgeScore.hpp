// no-networkit-format
/*
 * RandomEdgeScore.hpp
 *
 *  Created on: 11.08.2014
 *      Author: Gerd Lindner
 */

#ifndef NETWORKIT_SPARSIFICATION_RANDOM_EDGE_SCORE_HPP_
#define NETWORKIT_SPARSIFICATION_RANDOM_EDGE_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

/**
 * Generates a random edge attribute. Each edge is assigned a random value in [0,1].
 */
class RandomEdgeScore final : public EdgeScore<double> {

public:

    /**
     * Creates a new instance of the Random edge score.
     */
    RandomEdgeScore(const Graph& G);

    double score(edgeid eid) override;
    double score(node u, node v) override;
    void run() override;
};

}
/* namespace NetworKit */
#endif // NETWORKIT_SPARSIFICATION_RANDOM_EDGE_SCORE_HPP_
