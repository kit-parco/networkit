/*
 * OverlappingCommunityDetectionAlgorithm.hpp
 *
 *  Created on: 14.12.2020
 *      Author: John Gelhausen
 */
// networkit-format

#ifndef NETWORKIT_COMMUNITY_OVERLAPPING_COMMUNITY_DETECTION_ALGORITHM_HPP_
#define NETWORKIT_COMMUNITY_OVERLAPPING_COMMUNITY_DETECTION_ALGORITHM_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Abstract base class for overlapping community detection/graph clustering algorithms.
 */
class OverlappingCommunityDetectionAlgorithm : public Algorithm {
public:
    /**
     * An overlapping community detection algorithm operates on a graph, so the constructor expects
     * a graph.
     *
     * @param[in] G input graph
     */
    OverlappingCommunityDetectionAlgorithm(const Graph &G);

    /** Default destructor */
    ~OverlappingCommunityDetectionAlgorithm() override = default;

    /**
     * Apply algorithm to graph
     */
    void run() override = 0;

    /**
     * Returns the result of the run method or throws an error, if the algorithm hasn't run yet.
     * @return cover of the node set
     */
    const Cover &getCover() const;

protected:
    const Graph *G;
    Cover result;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_OVERLAPPING_COMMUNITY_DETECTION_ALGORITHM_HPP_
