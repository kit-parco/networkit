// no-networkit-format
/*
 * EdgeScoreLinearizer.hpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef NETWORKIT_EDGESCORES_EDGE_SCORE_LINEARIZER_HPP_
#define NETWORKIT_EDGESCORES_EDGE_SCORE_LINEARIZER_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

class EdgeScoreLinearizer final : public EdgeScore<double> {

private:
    const std::vector<double>* attribute;
    bool inverse;

public:
    EdgeScoreLinearizer(const Graph& graph, const std::vector<double>& attribute, bool inverse = false);

    double score(edgeid eid) override;
    double score(node u, node v) override;
    void run() override;

};

} // namespace NetworKit

#endif // NETWORKIT_EDGESCORES_EDGE_SCORE_LINEARIZER_HPP_
