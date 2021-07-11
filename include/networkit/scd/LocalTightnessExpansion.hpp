#ifndef NETWORKIT_SCD_LOCAL_TIGHTNESS_EXPANSION_HPP_
#define NETWORKIT_SCD_LOCAL_TIGHTNESS_EXPANSION_HPP_

#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {

/**
 * The Local Tightness Expansion (LTE) algorithm.
 *
 * The algorithm can handle weighted graphs.
 *
 * This is the local community expansion algorithm described in
 *
 * Huang, J., Sun, H., Liu, Y., Song, Q., & Weninger, T. (2011).
 * Towards Online Multiresolution Community Detection in Large-Scale Networks.
 * PLOS ONE, 6(8), e23829.
 * https://doi.org/10.1371/journal.pone.0023829
 */
class LocalTightnessExpansion : public SelectiveCommunityDetector {
private:
    double alpha;

public:
    /**
     * Constructs the Local Tightness Expansion algorithm.
     *
     * @param[in] G The graph to detect communities on
     * @param[in] alpha Tightness coefficient - smaller values lead to larger communities
     */
    LocalTightnessExpansion(const Graph &g, double alpha = 1.0);

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

#endif // NETWORKIT_SCD_LOCAL_TIGHTNESS_EXPANSION_HPP_
