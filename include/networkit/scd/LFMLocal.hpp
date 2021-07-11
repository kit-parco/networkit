#ifndef NETWORKIT_SCD_LFM_LOCAL_HPP_
#define NETWORKIT_SCD_LFM_LOCAL_HPP_

#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {

/**
 * Local version of the LFM algorithm
 *
 * This is the local community expansion as introduced in:
 *
 * Lancichinetti, A., Fortunato, S., & Kertész, J. (2009).
 * Detecting the overlapping and hierarchical community structure in complex networks.
 * New Journal of Physics, 11(3), 033015.
 * https://doi.org/10.1088/1367-2630/11/3/033015
 *
 * Their algorithm detects overlapping communities by repeatedly
 * executing this algorithm for a random seed node that has not yet
 * been assigned to any community.
 *
 * The algorithm has a resolution parameter alpha. A natural choice
 * for alpha is 1, the paper states that values below 0.5 usually
 * give a community containing the whole graph while values larger
 * than 2 recover the smallest communities.
 */
class LFMLocal : public SelectiveCommunityDetector {

public:
    /**
     * Construct the LFMLocal algorithm.
     *
     * @param G The graph to find a community on.
     * @param alpha The resolution parameter.
     */
    LFMLocal(const Graph &g, double alpha = 1.0);

    using SelectiveCommunityDetector::expandOneCommunity;

    /**
     * Expand a set of nodes into a single community.
     *
     * @param s The set of seed nodes.
     * @return The found community.
     */
    std::set<node> expandOneCommunity(const std::set<node> &s) override;

protected:
    const double alpha;
};

} // namespace NetworKit

#endif // NETWORKIT_SCD_LFM_LOCAL_HPP_
