// no-networkit-format
/*
 * EdgeScoreNormalizer.hpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef NETWORKIT_EDGESCORES_EDGE_SCORE_NORMALIZER_HPP_
#define NETWORKIT_EDGESCORES_EDGE_SCORE_NORMALIZER_HPP_

#include <networkit/edgescores/EdgeScore.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

template <typename A>
class EdgeScoreNormalizer final : public EdgeScore<double> {

public:
    EdgeScoreNormalizer(const Graph &G, const std::vector<A> &score, bool invert = false, double lower = 0, double upper = 1.0);

    double score(edgeid eid) override;
    double score(node u, node v) override;
    void run() override;

private:
    const std::vector<A> *input;
    bool invert;
    double lower, upper;
};

}

#endif // NETWORKIT_EDGESCORES_EDGE_SCORE_NORMALIZER_HPP_
