/* GCE.hpp
 *
 * Created on: 06.05.2013
 * Author: cls
 */

#ifndef NETWORKIT_SCD_GCE_HPP_
#define NETWORKIT_SCD_GCE_HPP_

#include <unordered_set>

#include <tlx/define/deprecated.hpp>
#include <networkit/auxiliary/SetIntersector.hpp>
#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {

/**
 * The Greedy Community Expansion algorithm.
 *
 * Greedily adds nodes from the shell to improve community quality.
 */
class GCE : public SelectiveCommunityDetector {

public:
    GCE(const Graph &g, std::string objective);

    /**
     * Expand a node into a community.
     *
     * Deprecated, use GCE::expandOneCommunity instead.
     *
     * @param seed Seed node for which a community is to be found.
     *
     * @return Set of nodes that makes up the best community found around the @a seed node.
     */
    std::set<node> TLX_DEPRECATED(expandSeed(node seed));

    /**
     * @param[in]  seeds  seed nodes
     *
     * @param[out]        community as a set of nodes
     */
    std::set<node> expandOneCommunity(const std::set<node> &seeds) override;

    // inherit method from parent class.
    using SelectiveCommunityDetector::expandOneCommunity;

private:
    std::string objective; // name of objective function
};

} /* namespace NetworKit */
#endif // NETWORKIT_SCD_GCE_HPP_
