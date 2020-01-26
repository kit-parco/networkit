/*
 * SingleCommunityCleanup.hpp
 *
 * Created: 2019-09-26
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_COMMUNITY_CLEANUP_SINGLE_COMMUNITY_CLEAN_UP_HPP_
#define NETWORKIT_COMMUNITY_CLEANUP_SINGLE_COMMUNITY_CLEAN_UP_HPP_

#include <networkit/auxiliary/SparseVector.hpp>
#include <networkit/community/cleanup/SignificanceCalculator.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Provides a way to clean up a detected community to improve its quality.
 * Based on the statistical significance of the communities proposed by Lancichinetti et al: Finding
 * statistically significant communities in networks.
 */
class SingleCommunityCleanUp {
public:
    using Community = std::vector<node>;

    explicit SingleCommunityCleanUp(const Graph &graph,
                                    const StochasticDistributionCalculator &stochasticDistribution,
                                    double scoreThreshold = 0.1,
                                    double significanceThreshold = 0.1,
                                    double minOverlapRatio = 0.5);

    /**
     * Clean a community using statistical significance. If the community is not considered
     * significant, an empty community is returned.
     * @param inputCommunity community to clean
     * @return cleaned community
     */
    Community clean(const Community &inputCommunity);

private:
    struct ScoreStruct {
        double rScore;
        node nodeId;

        explicit ScoreStruct(double score, node nodeId)
                : rScore(score), nodeId(nodeId) {};
    };

    const Graph &graph;
    SignificanceCalculator significanceCalculator;
    // threshold to decide if a node is significant
    const double significanceThreshold;
    // threshold to discard candidates because they will most likely not be significant
    const double scoreThreshold;
    // threshold to discard communities if they changed too much
    const double minOverlapRatio;

    std::vector<node> candidates;
    SparseVector<count> edgesToCommunity;
    SparseVector<int> isInCommunity;
    SparseVector<int> isInOriginalCommunity;
    SparseVector<int> isCandidate;
    count outgoingCommunityStubs;   // number of edges that have one endpoint inside the community and one outside
    count totalCommunityStubs;      // the sum of the degrees of the nodes in the community
    count externalNodes;
    count externalStubs;

    void getCandidatesAndSetUpCalculation(bool onlyUseOriginalCommunity);

    Community
    calculateSignificantNodes(const Community &inputCommunity, bool onlyUseOriginalCommunity);

    std::vector<ScoreStruct> calculateCandidateScores();

    std::vector<ScoreStruct> calculateInternalScores();

    std::vector<node> findSignificantCandidates(const std::vector<ScoreStruct> &scores);

    void removeWorstNode(std::vector<ScoreStruct> &internalScores);

    void reset();

    bool haveSmallOverlap(const Community &inputCommunity, const Community &cleanedCommunity);

    Community firstPhase(const Community &inputCommunity);

    Community secondPhase(const Community &firstPhaseResult);

    void removeNodeFromCommunity(node nodeToRemove);
};

} /* namespace NetworKit */

#endif // NETWORKIT_COMMUNITY_CLEANUP_SINGLE_COMMUNITY_CLEAN_UP_HPP_
