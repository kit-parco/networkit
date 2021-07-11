// no-networkit-format
/*
 * EdgeCut.hpp
 *
 *  Created on: Jun 20, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_COMMUNITY_EDGE_CUT_HPP_
#define NETWORKIT_COMMUNITY_EDGE_CUT_HPP_

#include <networkit/community/QualityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 */
class EdgeCut final : public QualityMeasure {
public:
    double getQuality(const Partition& zeta, const Graph& G) override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_EDGE_CUT_HPP_
