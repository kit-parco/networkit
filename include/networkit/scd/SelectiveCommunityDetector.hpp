/*
 * SelectiveCommunityDetector.hpp
 *
 *  Created on: 15.05.2013
 *      Author: cls
 */

#ifndef NETWORKIT_SCD_SELECTIVE_COMMUNITY_DETECTOR_HPP_
#define NETWORKIT_SCD_SELECTIVE_COMMUNITY_DETECTOR_HPP_

#include <map>
#include <set>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Base class for selective community detection algorithms, i.e., algorithms
 * that detect communities around a single node.
 */
class SelectiveCommunityDetector {

public:
    /**
     * Construct a selective community detector. This should only be called by child classes.
     *
     * @param G The graph for which communities shall be detected.
     */
    SelectiveCommunityDetector(const Graph &g);

    /**
     * Virtual default destructor to allow safe destruction of child classes.
     */
    virtual ~SelectiveCommunityDetector() = default;

    /**
     * Detect one community for each of the given seed nodes.
     *
     * The default implementation calls expandOneCommunity() for each of the seeds.
     *
     * @param seeds The list of seeds for which communities shall be detected.
     * @return a mapping from seed node to community (as a set of nodes)
     */
    virtual std::map<node, std::set<node>> run(const std::set<node> &seeds);

    /**
     * Detect a community for the given seed node.
     *
     * The default implementation calls expandOneCommunity(const
     * std::set<node>&) with a set of one node.
     *
     * @param seed The seed to find the community for.
     * @return The found community as set of node.
     */
    virtual std::set<node> expandOneCommunity(node seed);

    /**
     * Detect a single community for the given seed nodes.
     *
     * This is useful if you know multiple nodes that should be part of the
     * community. This method may throw an exception if the particular algorithm
     * does not support multiple seeds.
     *
     * @param seeds The seeds for the community.
     * @return The found community as set of nodes.
     */
    virtual std::set<node> expandOneCommunity(const std::set<node> &seeds) = 0;

protected:
    const Graph *g;
};

} /* namespace NetworKit */
#endif // NETWORKIT_SCD_SELECTIVE_COMMUNITY_DETECTOR_HPP_
