// no-networkit-format
/*
 * SimmelianOverlapScore.hpp
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#ifndef NETWORKIT_SPARSIFICATION_SIMMELIAN_OVERLAP_SCORE_HPP_
#define NETWORKIT_SPARSIFICATION_SIMMELIAN_OVERLAP_SCORE_HPP_

#include <networkit/sparsification/SimmelianScore.hpp>
#include <set>

namespace NetworKit {

/**
 * Calculates the Simmelian backbone (paramaetric variant) for a given input graph.
 */
class SimmelianOverlapScore final : public SimmelianScore {

public:

    /**
     * Creates a new instance of the parametric variant of the Simmelian Backbone calculator.
     * @param maxRank     the maximum rank that is considered for overlap calculation
     */
    SimmelianOverlapScore(const Graph& graph, const std::vector<count>& triangles, count maxRank);
    void run() override;

private:
    count maxRank;
};

}
/* namespace NetworKit */
#endif // NETWORKIT_SPARSIFICATION_SIMMELIAN_OVERLAP_SCORE_HPP_
