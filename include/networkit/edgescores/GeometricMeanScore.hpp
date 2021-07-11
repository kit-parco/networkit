// no-networkit-format
/*
 * GeometricMeanScore.hpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef NETWORKIT_EDGESCORES_GEOMETRIC_MEAN_SCORE_HPP_
#define NETWORKIT_EDGESCORES_GEOMETRIC_MEAN_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

class GeometricMeanScore final : public EdgeScore<double> {

private:
    const std::vector<double>* attribute;

public:
    GeometricMeanScore(const Graph& G, const std::vector<double>& attribute);
    double score(edgeid eid) override;
    double score(node u, node v) override;
    void run() override;
};

} // namespace NetworKit

#endif // NETWORKIT_EDGESCORES_GEOMETRIC_MEAN_SCORE_HPP_
