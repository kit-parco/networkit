// no-networkit-format
/*
 * LPDegreeOrdered.hpp
 *
 *  Created on: 24.09.2013
 *      Author: cls
 */

#ifndef NETWORKIT_COMMUNITY_LP_DEGREE_ORDERED_HPP_
#define NETWORKIT_COMMUNITY_LP_DEGREE_ORDERED_HPP_

#include <networkit/community/CommunityDetectionAlgorithm.hpp>

namespace NetworKit {

using label = index; // a label is the same as a cluster id

/**
 * @ingroup community
 * Label propagation-based community detection algorithm which
 * processes nodes in increasing order of node degree.
 */
class LPDegreeOrdered final : public CommunityDetectionAlgorithm {
private:
    count nIterations = 0; //!< number of iterations in last run

public:
    /**
     * Constructor to the degree ordered label propagation community detection algorithm.
     *
     * @param[in] G input graph
     */
    LPDegreeOrdered(const Graph& G);

    /**
     * Detect communities.
     */
    void run() override;

    /**
    * Get number of iterations in last run.
    *
    * @return Number of iterations.
    */
    count numberOfIterations();

    std::string TLX_DEPRECATED(toString() const override);

};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_LP_DEGREE_ORDERED_HPP_
