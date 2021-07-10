/*
 * LFM.hpp
 *
 *  Created on 10.12.2020
 *      Author: John Gelhausen
 */
// networkit-format

#ifndef NETWORKIT_COMMUNITY_LFM_HPP_
#define NETWORKIT_COMMUNITY_LFM_HPP_

#include <networkit/community/OverlappingCommunityDetectionAlgorithm.hpp>
#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {
/**
 * @ingroup community
 * This is the local community expansion algorithm as introduced in:
 *
 * Lancichinetti, A., Fortunato, S., & Kertész, J. (2009).
 * Detecting the overlapping and hierarchical community structure in complex networks.
 * New Journal of Physics, 11(3), 033015.
 * https://doi.org/10.1088/1367-2630/11/3/033015
 *
 * The LFM algorithm detects overlapping communities by repeatedly
 * executing a given SelectiveCommunityDetector algorithm
 * for different random seed nodes which have not yet been assigned to any community.
 * While this implementation allows the use of any local community detection algorithm, the behavior
 * of the original algorithm can be achieved using the LFMLocal local community detection algorithm.
 */
class LFM final : public OverlappingCommunityDetectionAlgorithm {
public:
    /**
     * @param G Input graph
     * @param scd The algorithm that is used to expand the random seed nodes to communities
     *
     */
    LFM(const Graph &G, SelectiveCommunityDetector &scd);

    /**
     * Detect communities
     */
    void run() override;

protected:
    SelectiveCommunityDetector *scd;
};
} /* namespace NetworKit */

#endif // NETWORKIT_COMMUNITY_LFM_HPP_
