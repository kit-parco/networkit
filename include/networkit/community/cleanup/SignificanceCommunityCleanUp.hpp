/*
 * SignificanceCommunityCleanUp.hpp
 *
 * Created: 2019-09-11
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_COMMUNITY_CLEANUP_SIGNIFICANCE_COMMUNITY_CLEAN_UP_HPP_
#define NETWORKIT_COMMUNITY_CLEANUP_SIGNIFICANCE_COMMUNITY_CLEAN_UP_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/community/cleanup/SignificanceCalculator.hpp>
#include <networkit/community/cleanup/SingleCommunityCleanUp.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * An algorithm that aims to improve the quality of (overlapping) communities, e.g. to clean up
 * the communities detected by a (overlapping) community detection algorithm.
 * Based on the statistical significance of the communities proposed by Lancichinetti et al: Finding
 * statistically significant communities in networks.
 */
class SignificanceCommunityCleanUp : public Algorithm {
public:
    /**
     * Constructor of the algorithm.
     * @param[in] graph	input graph
     * @param[inout] communities input communities that will be cleaned up
     * @param[in] distribution The stochastic distribution object to use for stochastic calculations
     */
    SignificanceCommunityCleanUp(const Graph &graph,
                                 std::vector<std::vector<node>> &communities,
                                 StochasticDistributionCalculator &distribution,
                                 double significanceThreshold = 0.1,
                                 double scoreThreshold = 0.1,
                                 double minOverlapRatio = 0.5,
                                 bool mergeDiscarded = true);

    void run() override;

    std::string toString() const override;

    bool isParallel() const override;

private:

    const Graph &graph;
    std::vector<std::vector<node>> &communities;
    std::vector<std::vector<node>> discardedCommunities;
    double significanceThreshold;
    double scoreThreshold;
    double minOverlapRatio;
    const bool mergeDiscarded;
    count maxCommunitySize;

    StochasticDistributionCalculator &stochasticDistribution;

    void cleanAllCommunities();

    void mergeDiscardedCommunities();

};
} /* namespace NetworKit */

#endif // NETWORKIT_COMMUNITY_CLEANUP_SIGNIFICANCE_COMMUNITY_CLEAN_UP_HPP_
