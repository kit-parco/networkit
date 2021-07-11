// no-networkit-format
/*
 * CommunityDetectionAlgorithm.hpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_COMMUNITY_COMMUNITY_DETECTION_ALGORITHM_HPP_
#define NETWORKIT_COMMUNITY_COMMUNITY_DETECTION_ALGORITHM_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Abstract base class for community detection/graph clustering algorithms.
 */
class CommunityDetectionAlgorithm : public Algorithm {
public:
    /**
     * A community detection algorithm operates on a graph, so the constructor expects a graph.
     *
     * @param[in] G input graph
     */
    CommunityDetectionAlgorithm(const Graph& G);

    /**
     * A community detection algorithm operates on a graph, so the constructor expects a graph.
     *
     * @param[in] G input graph
     * @param[in] baseClustering optional; the algorithm will start from the given clustering.
     */
    CommunityDetectionAlgorithm(const Graph &G, Partition baseClustering);

    /** Default destructor */
    ~CommunityDetectionAlgorithm() override = default;

    /**
     * Apply algorithm to graph
     */
    void run() override = 0;

    /**
     * Returns the result of the run method or throws an error, if the algorithm hasn't run yet.
     * @return partition of the node set
     */
    virtual Partition getPartition();

    /**
     * @return string representation of algorithm and parameters.
     */
    std::string TLX_DEPRECATED(toString() const override);

protected:
    const Graph* G;
    Partition result;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_COMMUNITY_DETECTION_ALGORITHM_HPP_
