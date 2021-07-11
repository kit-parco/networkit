// no-networkit-format
/*
 * Coverage.hpp
 *
 *  Created on: 02.02.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_COMMUNITY_COVERAGE_HPP_
#define NETWORKIT_COMMUNITY_COVERAGE_HPP_

#include <networkit/community/QualityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Coverage is the fraction of intra-cluster edges.
 */
class Coverage final: public QualityMeasure {
public:

    double getQuality(const Partition& zeta, const Graph& G) override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_COVERAGE_HPP_
