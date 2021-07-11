/*
 * ChangeCorrectedTriangleScore.hpp
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#ifndef NETWORKIT_SPARSIFICATION_CHANCE_CORRECTED_TRIANGLE_SCORE_HPP_
#define NETWORKIT_SPARSIFICATION_CHANCE_CORRECTED_TRIANGLE_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

class ChanceCorrectedTriangleScore final : public EdgeScore<double> {

public:
    ChanceCorrectedTriangleScore(const Graph &graph, const std::vector<count> &triangles);
    double score(edgeid eid) override;
    double score(node u, node v) override;
    void run() override;

private:
    const std::vector<count> *triangles;
};

} // namespace NetworKit
#endif // NETWORKIT_SPARSIFICATION_CHANCE_CORRECTED_TRIANGLE_SCORE_HPP_
