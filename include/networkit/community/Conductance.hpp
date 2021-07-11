// no-networkit-format
/*
 * Conductance.hpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#ifndef NETWORKIT_COMMUNITY_CONDUCTANCE_HPP_
#define NETWORKIT_COMMUNITY_CONDUCTANCE_HPP_

#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/community/EdgeCut.hpp>
#include <networkit/community/QualityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Compute conductance of a 2-partition, i.e. cut size over volume of smaller set (smaller in
 * terms of volume).
 */
class Conductance final : public QualityMeasure {
public:
    /**
     * @return Conductance of 2-partition @a zeta in graph @a G.
     * Requires cluster IDs to be either 0 or 1.
     */
    double getQuality(const Partition& zeta, const Graph& G) override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_CONDUCTANCE_HPP_
