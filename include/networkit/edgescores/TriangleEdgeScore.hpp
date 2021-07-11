/*
 * TriangleEdgeScore.hpp
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner, Michael Hamann
 */

#ifndef NETWORKIT_EDGESCORES_TRIANGLE_EDGE_SCORE_HPP_
#define NETWORKIT_EDGESCORES_TRIANGLE_EDGE_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

/**
 * A parallel triangle counting implementation based on ideas in [0].
 *
 * This implementation avoids less work but needs therefore also less checks and
 * is (apart from a fast initialization) parallel without any locks. In experiments
 * this implementation seems to be fast both in non-parallel as well as in parallel
 * settings. With only one thread its performance is similar to the sequential
 * ChibaNishizekiTriangleEdgeScore in NetworKit.
 *
 * [0] Triangle Listing Algorithms: Back from the Diversion
 * Mark Ortmann and Ulrik Brandes * 2014 Proceedings of the Sixteenth Workshop on Algorithm
 * Engineering and Experiments (ALENEX). 2014, 1-8
 */
class TriangleEdgeScore final : public EdgeScore<count> {

public:
    TriangleEdgeScore(const Graph &G);
    count score(edgeid eid) override;
    count score(node u, node v) override;
    void run() override;
};

} /* namespace NetworKit */

#endif // NETWORKIT_EDGESCORES_TRIANGLE_EDGE_SCORE_HPP_
