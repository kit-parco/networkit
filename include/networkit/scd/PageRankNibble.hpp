/*
 * PageRankNibble.hpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#ifndef NETWORKIT_SCD_PAGE_RANK_NIBBLE_HPP_
#define NETWORKIT_SCD_PAGE_RANK_NIBBLE_HPP_

#include <set>

#include <tlx/define/deprecated.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {

/**
 * Variant of PageRank-Nibble algorithm due to Andersen, Chung and Lang.
 * Paper: Local Graph Partitioning using PageRank Vectors.
 * URL: http://www.math.ucsd.edu/~fan/wp/localpartition.pdf
 * Simplifications according to D. Gleich's code at URL https://gist.github.com/dgleich/6201856.
 */
class PageRankNibble final : public SelectiveCommunityDetector {

    double alpha;
    double epsilon;

    std::set<node> bestSweepSet(std::vector<std::pair<node, double>> &pr);

public:
    /**
     * @param Graph @a g for which PageRank-Nibble is to be performed. Is treated as
     * unweighted in current implementation.
     * @param alpha Loop probability of random walk; smaller values tend to produce larger
     *        communities.
     * @param eps Tolerance threshold for approximation of PageRank vectors.
     */
    PageRankNibble(const Graph &g, double alpha, double epsilon);

    ~PageRankNibble() override = default;

    /**
     * Expand a node into a community.
     *
     * Deprecated, use PageRankNibble::expandOneCommunity instead.
     *
     * @param seed Seed node for which a community is to be found.
     *
     * @return Set of nodes that makes up the best community found around the @a seed node.
     */
    std::set<node> TLX_DEPRECATED(expandSeed(node seed));

    /**
     * @param seeds Seed nodes for which a community is to be found.
     *
     * @return Set of nodes that makes up the best community found around the @a seeds nodes.
     */
    std::set<node> expandOneCommunity(const std::set<node> &seeds) override;

    // inherit method from parent class.
    using SelectiveCommunityDetector::expandOneCommunity;
};

} /* namespace NetworKit */
#endif // NETWORKIT_SCD_PAGE_RANK_NIBBLE_HPP_
