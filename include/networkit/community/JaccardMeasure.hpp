// no-networkit-format
/*
 * JaccardMeasure.hpp
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_COMMUNITY_JACCARD_MEASURE_HPP_
#define NETWORKIT_COMMUNITY_JACCARD_MEASURE_HPP_

#include <networkit/community/DissimilarityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 */
class JaccardMeasure final: public DissimilarityMeasure {

public:

    double getDissimilarity(const Graph &G, const Partition &zeta, const Partition &eta) override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_JACCARD_MEASURE_HPP_
