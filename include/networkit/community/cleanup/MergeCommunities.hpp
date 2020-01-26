/*
 * MergeCommunities.hpp
 *
 * Created: 2019-09-26
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_COMMUNITY_CLEANUP_MERGE_COMMUNITIES_HPP_
#define NETWORKIT_COMMUNITY_CLEANUP_MERGE_COMMUNITIES_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/community/cleanup/SingleCommunityCleanUp.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * Merge not significant communities to find significant communities, using the statistical
 * significance.
 */
class MergeCommunities : public Algorithm {
public:
    MergeCommunities(const Graph &graph,
                     std::vector<std::vector<node>> discardedCommunities,
                     StochasticDistributionCalculator &stochasticDistribution,
                     double significanceThreshold = 0.1,
                     double scoreThreshold = 0.1,
                     double minOverlapRatio = 0.5,
                     count maxCommunitySize = none);

    void run() override;

    const std::vector<std::vector<node>>& getCleanedCommunities();

    std::string toString() const override;

    bool isParallel() const override;

private:
    const Graph &graph;
    std::vector<std::vector<node>> discardedCommunities;
    StochasticDistributionCalculator &stochasticDistribution;
    SignificanceCalculator significanceCalculator;
    double significanceThreshold;
    double scoreThreshold;
    double minOverlapRatio;
    std::vector<std::vector<node>> cleanedCommunities;
    Graph discardedCommunitiesGraph;
    Partition mergedCommunities;
    std::vector<count> outgoingGroupStubs;
    std::vector<count> totalGroupStubs;
    count totalStubs;
    const count maxCommunitySize;

    void createDiscardedCommunitiesGraph();

    void tryToMergeCommunities();

    void checkMergedCommunities();

    bool tryLocalMove(node u, SparseVector<edgeweight> &neighborWeights);
};


} /* namespace NetworKit */


#endif // NETWORKIT_COMMUNITY_CLEANUP_MERGE_COMMUNITIES_HPP_
