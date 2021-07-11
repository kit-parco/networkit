// no-networkit-format
/*
 * EdgeScoreBlender.hpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef NETWORKIT_EDGESCORES_EDGE_SCORE_BLENDER_HPP_
#define NETWORKIT_EDGESCORES_EDGE_SCORE_BLENDER_HPP_

#include <networkit/edgescores/EdgeScore.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class EdgeScoreBlender final : public EdgeScore<double> {

public:

    EdgeScoreBlender(const Graph &G, const std::vector<double> &attribute0, const std::vector<double> &attribute1, const std::vector<bool> &selection);

    double score(edgeid eid) override;
    double score(node u, node v) override;
    void run() override;

private:
    const std::vector<double> *attribute0, *attribute1;
    const std::vector<bool> *selection;
};

} // namespace NetworKit

#endif // NETWORKIT_EDGESCORES_EDGE_SCORE_BLENDER_HPP_
