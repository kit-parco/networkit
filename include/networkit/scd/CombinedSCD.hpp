#ifndef NETWORKIT_SCD_COMBINED_SCD_HPP_
#define NETWORKIT_SCD_COMBINED_SCD_HPP_

#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {
/**
 * Helper for combining two selective community detection algorithms.
 *
 * This allows to combine two SCD algorithms such that first the first
 * algorithm is executed on the given seed(s) and then the result of
 * the first algorithm is used as the seed set of the second algorithm.
 * This is particularly useful for seeding an algorithm with a clique.
 *
 * @author Michael Hamann <michael.hamann@kit.edu>
 */
class CombinedSCD : public SelectiveCommunityDetector {
public:
    /**
     * Initialize the combined algorithm with the given graph and the given two algorithms.
     *
     * @param G The graph to work on.
     * @param first The first algorithm that is run with the given seed(s).
     * @param secodn The second algorithm that is run with the result of the first algorithm.
     */
    CombinedSCD(const Graph &g, SelectiveCommunityDetector &first,
                SelectiveCommunityDetector &second);

    /**
     * Expand a community with the given seed node using first and second algorithm.
     *
     * @param s The seed node to start with.
     * @return The community found by the second algorithm.
     */
    std::set<node> expandOneCommunity(node s) override;

    /**
     * Expand a community with the given seed nodes using first and second algorithm.
     *
     * @param s The seed nodes to start with.
     * @return The community found by the second algorithm.
     */
    std::set<node> expandOneCommunity(const std::set<node> &s) override;

protected:
    SelectiveCommunityDetector &first;
    SelectiveCommunityDetector &second;
};
} // namespace NetworKit

#endif // NETWORKIT_SCD_COMBINED_SCD_HPP_
