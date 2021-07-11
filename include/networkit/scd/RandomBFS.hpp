#ifndef NETWORKIT_SCD_RANDOM_BFS_HPP_
#define NETWORKIT_SCD_RANDOM_BFS_HPP_

#include <networkit/scd/SelectiveCommunityDetector.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {

/*
 * The random BFS community detection baseline:
 * finds a community around a seed node with a given size using a prefix of a
 * BFS order that is selected at random among all possible prefixes.
 */
class RandomBFS : public SelectiveCommunityDetector {

public:
    RandomBFS(const Graph &g, const Cover &c);

    /**
     * @param[in]	s	seed node
     *
     * @param[out]		community as a set of nodes
     */
    std::set<node> expandOneCommunity(const std::set<node> &s) override;
    using SelectiveCommunityDetector::expandOneCommunity;

protected:
    const Cover *c; // ground truth communities
    std::map<index, count> subsetSizes;
};

} /* namespace NetworKit */
#endif // NETWORKIT_SCD_RANDOM_BFS_HPP_
