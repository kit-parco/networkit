#ifndef NETWORKIT_SCD_TCE_HPP_
#define NETWORKIT_SCD_TCE_HPP_

#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {

/**
 * Triangle-based community expansion
 *
 * The algorithm can handle weighted graphs.
 */
class TCE : public SelectiveCommunityDetector {
private:
    bool refine;
    bool useJaccard;

public:
    /**
     * Construct a TCE object.
     *
     * @param[in] G The graph to detect communities on
     * @param[in] refine If nodes shall be removed again if this makes the community better
     * @param[in] useJaccard use jaccard index for weight calculation.
     */
    TCE(const Graph &g, bool refine = false, bool useJaccard = false);

    /**
     * Expands a set of seed nodes into a community.
     *
     * @param[in] s The seed nodes
     * @return A community of the seed nodes
     */
    std::set<node> expandOneCommunity(const std::set<node> &s) override;

    using SelectiveCommunityDetector::expandOneCommunity;
};

} // namespace NetworKit

#endif // NETWORKIT_SCD_TCE_HPP_
