// no-networkit-format
/*
 * ForestFireScore.hpp
 *
 *  Created on: 26.08.2014
 *      Author: Gerd Lindner
 */

#ifndef NETWORKIT_SPARSIFICATION_FOREST_FIRE_SCORE_HPP_
#define NETWORKIT_SPARSIFICATION_FOREST_FIRE_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

/**
 * Based on the Forest Fire algorithm introduced by Leskovec et al.
 * The burn frequency of the edges is used as edge score.
 */
class ForestFireScore final : public EdgeScore<double> {

public:

    ForestFireScore(const Graph& graph, double pf, double targetBurntRatio);
    double score(edgeid eid) override;
    double score(node u, node v) override;
    void run() override;

private:
    double pf;
    double targetBurntRatio;

};

}
/* namespace NetworKit */
#endif // NETWORKIT_SPARSIFICATION_FOREST_FIRE_SCORE_HPP_
