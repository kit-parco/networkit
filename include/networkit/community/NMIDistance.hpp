// no-networkit-format
/*
 * NMIDistance.hpp
 *
 *  Created on: 30.04.2013
 *      Author: cls
 */

#ifndef NETWORKIT_COMMUNITY_NMI_DISTANCE_HPP_
#define NETWORKIT_COMMUNITY_NMI_DISTANCE_HPP_

#include <networkit/community/DissimilarityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * NMIDistance quantifies the dissimilarity between two clusterings using
 * Normalized Mutual Information.
 *
 */
class NMIDistance final: public DissimilarityMeasure {

public:

    double getDissimilarity(const Graph& G, const Partition& zeta, const Partition& eta) override;

};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_NMI_DISTANCE_HPP_
