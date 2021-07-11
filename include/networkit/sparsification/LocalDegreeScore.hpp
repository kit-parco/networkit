// no-networkit-format
/*
 * LocalDegreeScore.hpp
 *
 *  Created on: 28.08.2014
 *      Author: Gerd Lindner
 */

#ifndef NETWORKIT_SPARSIFICATION_LOCAL_DEGREE_SCORE_HPP_
#define NETWORKIT_SPARSIFICATION_LOCAL_DEGREE_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

/**
 * Local Degree sparsification method.
 * See 'Structure-Preserving Sparsification of Social Networks' by Lindner, Staudt, Hamann.
 */
class LocalDegreeScore final : public EdgeScore<double> {

public:

    LocalDegreeScore(const Graph& G);
    void run() override;
    double score(edgeid eid) override;
    double score(node u, node v) override;

};

}
/* namespace NetworKit */
#endif // NETWORKIT_SPARSIFICATION_LOCAL_DEGREE_SCORE_HPP_
