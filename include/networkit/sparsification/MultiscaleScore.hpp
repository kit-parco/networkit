// no-networkit-format
/*
 * MultiscaleScore.hpp
 *
 *  Created on: 20.06.2014
 *      Author: Gerd Lindner
 */

#ifndef NETWORKIT_SPARSIFICATION_MULTISCALE_SCORE_HPP_
#define NETWORKIT_SPARSIFICATION_MULTISCALE_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

/**
 * Calculates the multiscale edge score for a given graph. Each edge is
 * assigned the maximum filter value in [0,1] for which the edge will be contained
 * in the multiscale backbone.
 *
 * See "Extracting the multiscale backbone of complex weighted networks" by Serrano et al.
 */
class MultiscaleScore final : public EdgeScore<double> {

public:

    MultiscaleScore(const Graph& graph, const std::vector<double>& attribute);
    double score(edgeid eid) override;
    double score(node u, node v) override;
    void run() override;
    double getProbability(count degree, edgeweight normalizedWeight);

private:
    const std::vector<double>* attribute;
};

}
/* namespace NetworKit */

#endif // NETWORKIT_SPARSIFICATION_MULTISCALE_SCORE_HPP_
