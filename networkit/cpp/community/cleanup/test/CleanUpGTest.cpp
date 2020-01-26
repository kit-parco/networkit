/*
 * CleanUpGTest.cpp
 *
 * Created: 2019-09-18
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>
#include <cmath>

#include <networkit/community/EgoSplitting.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/community/cleanup/SignificanceCommunityCleanUp.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/community/cleanup/MergeCommunities.hpp>

namespace NetworKit {

class CleanupGTest : public testing::Test {
public:
    void SetUp() {
        Aux::Random::setSeed(435913, false);
    }
};

TEST_F(CleanupGTest, testSingleCommunityCleanUp) {
    METISGraphReader graphReader;
    Graph G = graphReader.read("input/erdos_renyi_200_0.05.graph");
    // Create clique
    count cliqueSize = 6;
    for (node u = 0; u < cliqueSize; ++u) {
        for (node v = u + 1; v < cliqueSize; ++v) {
            if (!G.hasEdge(u, v))
                G.addEdge(u, v, defaultEdgeWeight);
        }
    }
    std::vector<node> expectedCommunity;
    for (node u = 0; u < cliqueSize; ++u)
        expectedCommunity.push_back(u);
    // Create a community for the clique, but exclude a node and include weakly connected ones
    std::vector<node> testCommunity;
    count excludeCliqueMembers = 1;
    count addWeaklyConnected = 3;
    for (node u = excludeCliqueMembers; u < cliqueSize + addWeaklyConnected; ++u)
        testCommunity.push_back(u);
    StochasticDistributionCalculator dist(2 * G.numberOfEdges() + G.numberOfNodes());
    SingleCommunityCleanUp singleCommunityCleanUp(G, dist);

    std::vector<node> cleanedCommunity = singleCommunityCleanUp.clean(testCommunity);
    std::sort(cleanedCommunity.begin(), cleanedCommunity.end());

    EXPECT_EQ(cleanedCommunity, expectedCommunity);
}

TEST_F(CleanupGTest, testMergeDiscarded) {
    METISGraphReader graphReader;
    Graph G = graphReader.read("input/erdos_renyi_200_0.05.graph");
    // Create clique
    count cliqueSize = 8;
    for (node u = 0; u < cliqueSize; ++u) {
        for (node v = u + 1; v < cliqueSize; ++v) {
            if (!G.hasEdge(u, v))
                G.addEdge(u, v, defaultEdgeWeight);
        }
    }
    std::vector<node> expectedCommunity;
    for (node u = 0; u < cliqueSize; ++u)
        expectedCommunity.push_back(u);
    std::vector<std::vector<node>> discardedCommunitites;
    // Break clique into 4 discarded communities
    discardedCommunitites.push_back({0, 1});
    discardedCommunitites.push_back({2, 3});
    discardedCommunitites.push_back({4, 5});
    discardedCommunitites.push_back({6, 7});
    // Add some bad communities
    discardedCommunitites.push_back({10, 11, 12, 13});
    discardedCommunitites.push_back({15, 16});
    discardedCommunitites.push_back({18});
    discardedCommunitites.push_back({19});
    StochasticDistributionCalculator dist(2 * G.numberOfEdges() + G.numberOfNodes());
    MergeCommunities mergeCommunities(G, discardedCommunitites, dist);

    mergeCommunities.run();
    auto cleanedCommunities = mergeCommunities.getCleanedCommunities();

    EXPECT_EQ(cleanedCommunities.size(), 1);
    auto cleanedCommunity = cleanedCommunities.front();
    std::sort(cleanedCommunity.begin(), cleanedCommunity.end());
    EXPECT_EQ(cleanedCommunity, expectedCommunity);
}

TEST_F(CleanupGTest, testCleanUp) {
    METISGraphReader graphReader;
    Graph G = graphReader.read("input/10_clusters.graph");
    node isolatedNode = G.addNode();
    std::map<std::string, std::string> parameters;
    parameters["Extend EgoNet Strategy"] = "None";
    EgoSplitting algo(G);
    algo.setParameters(parameters);
    algo.run();
    Cover cover = algo.getCover();
    // Add bad communities
    cover.addSubset({1});
    cover.addSubset({2, isolatedNode});
    std::vector<std::vector<node>> communities(cover.upperBound());
    cover.forEntries([&](node u, const std::set<index>& coms) {
        for (index s : coms) {
            communities[s].push_back(u);
        }
    });
    count numberOfInputCommunities = communities.size();

    StochasticDistributionCalculator dist(2 * G.numberOfEdges() + G.numberOfNodes());
    SignificanceCommunityCleanUp cleanUp(G, communities, dist, 0.1, 0.1, 0.5);
    cleanUp.run();

    EXPECT_TRUE(communities.size() <= numberOfInputCommunities);
    // Communities of size 1 should be discarded
    for (const auto& c : communities) {
        EXPECT_GT(c.size(), 1);
    }
    count notEmptyComms = 0;
    for (const auto& c : communities) {
        notEmptyComms += (c.size() > 1);
    }
    EXPECT_GE(notEmptyComms, 6);
    std::vector<index> badComm = {2, isolatedNode};
    for (auto &comm : communities) {
        std::sort(comm.begin(), comm.end());
        EXPECT_NE(comm, badComm);
    }
}

TEST_F(CleanupGTest, testBinomialCoeff) {
    StochasticDistributionCalculator stoch(10);

    EXPECT_NEAR(1, stoch.binomialCoefficient(7, 0), 1e-6);
    EXPECT_NEAR(7, stoch.binomialCoefficient(7, 1), 1e-6);
    EXPECT_NEAR(21, stoch.binomialCoefficient(7, 2), 1e-6);
    EXPECT_NEAR(35, stoch.binomialCoefficient(7, 3), 1e-6);
    EXPECT_NEAR(35, stoch.binomialCoefficient(7, 4), 1e-6);
    EXPECT_NEAR(21, stoch.binomialCoefficient(7, 5), 1e-6);
    EXPECT_NEAR(7, stoch.binomialCoefficient(7, 6), 1e-6);
    EXPECT_NEAR(1, stoch.binomialCoefficient(7, 7), 1e-6);
}

TEST_F(CleanupGTest, testBinomialDist) {
    StochasticDistributionCalculator stoch(10);

    EXPECT_NEAR(stoch.binomialDistribution(0.3, 7, 0), 0.0823543, 1e-7);
    EXPECT_NEAR(stoch.binomialDistribution(0.3, 7, 1), 0.247063, 1e-6);
    EXPECT_NEAR(stoch.binomialDistribution(0.3, 7, 2), 0.317652, 1e-6);
    EXPECT_NEAR(stoch.binomialDistribution(0.3, 7, 3), 0.226894, 1e-6);
    EXPECT_NEAR(stoch.binomialDistribution(0.3, 7, 4), 0.0972405, 1e-7);
    EXPECT_NEAR(stoch.binomialDistribution(0.3, 7, 5), 0.0250047, 1e-7);
    EXPECT_NEAR(stoch.binomialDistribution(0.3, 7, 6), 0.0035721, 1e-8);
    EXPECT_NEAR(stoch.binomialDistribution(0.3, 7, 7), 0.0002187, 1e-9);
}

TEST_F(CleanupGTest, testRightCumBinom) {
    StochasticDistributionCalculator stoch(10);

    EXPECT_NEAR(stoch.rightCumulativeBinomial(0.3, 5, 0), 1, 1e-6);
    EXPECT_NEAR(stoch.rightCumulativeBinomial(0.3, 5, 1), 0.83193, 1e-6);
    EXPECT_NEAR(stoch.rightCumulativeBinomial(0.3, 5, 2), 0.47178, 1e-6);
    EXPECT_NEAR(stoch.rightCumulativeBinomial(0.3, 5, 3), 0.16308, 1e-6);
    EXPECT_NEAR(stoch.rightCumulativeBinomial(0.3, 5, 4), 0.03078, 1e-7);
    EXPECT_NEAR(stoch.rightCumulativeBinomial(0.3, 5, 5), 0.00243, 1e-8);
}

TEST_F(CleanupGTest, testLeftCumBinom) {
    StochasticDistributionCalculator stoch(10);

    EXPECT_NEAR(stoch.leftCumulativeBinomial(0.3, 5, 0), 0.16807, 1e-6);
    EXPECT_NEAR(stoch.leftCumulativeBinomial(0.3, 5, 1), 0.52822, 1e-6);
    EXPECT_NEAR(stoch.leftCumulativeBinomial(0.3, 5, 2), 0.83692, 1e-6);
    EXPECT_NEAR(stoch.leftCumulativeBinomial(0.3, 5, 3), 0.96922, 1e-6);
    EXPECT_NEAR(stoch.leftCumulativeBinomial(0.3, 5, 4), 0.99757, 1e-6);
    EXPECT_NEAR(stoch.leftCumulativeBinomial(0.3, 5, 5), 1, 1e-6);
}

TEST_F(CleanupGTest, testHyperDist) {
    StochasticDistributionCalculator stoch(10);

    EXPECT_NEAR(stoch.hypergeometricDistribution(10, 5, 5, 0), 0.003968254, 1e-8);
    EXPECT_NEAR(stoch.hypergeometricDistribution(10, 5, 5, 1), 0.099206349, 1e-7);
    EXPECT_NEAR(stoch.hypergeometricDistribution(10, 5, 5, 2), 0.396825397, 1e-6);
    EXPECT_NEAR(stoch.hypergeometricDistribution(10, 5, 5, 3), 0.396825397, 1e-6);
    EXPECT_NEAR(stoch.hypergeometricDistribution(10, 5, 5, 4), 0.099206349, 1e-7);
    EXPECT_NEAR(stoch.hypergeometricDistribution(10, 5, 5, 5), 0.003968254, 1e-8);
}

TEST_F(CleanupGTest, testRightCumHyper) {
    StochasticDistributionCalculator stoch(10);

    EXPECT_NEAR(stoch.rightCumulativeHypergeometric(10, 5, 5, 0), 1, 1e-6);
    EXPECT_NEAR(stoch.rightCumulativeHypergeometric(10, 5, 5, 1), 0.996031746, 1e-6);
    EXPECT_NEAR(stoch.rightCumulativeHypergeometric(10, 5, 5, 2), 0.896825397, 1e-6);
    EXPECT_NEAR(stoch.rightCumulativeHypergeometric(10, 5, 5, 3), 0.5, 1e-6);
    EXPECT_NEAR(stoch.rightCumulativeHypergeometric(10, 5, 5, 4), 0.103174603, 1e-6);
    EXPECT_NEAR(stoch.rightCumulativeHypergeometric(10, 5, 5, 5), 0.003968254, 1e-8);
}

TEST_F(CleanupGTest, testLeftCumHyper) {
    StochasticDistributionCalculator stoch(10);

    EXPECT_NEAR(stoch.leftCumulativeHypergeometric(10, 5, 5, 0), 0.003968254, 1e-8);
    EXPECT_NEAR(stoch.leftCumulativeHypergeometric(10, 5, 5, 1), 0.103174603, 1e-6);
    EXPECT_NEAR(stoch.leftCumulativeHypergeometric(10, 5, 5, 2), 0.5, 1e-6);
    EXPECT_NEAR(stoch.leftCumulativeHypergeometric(10, 5, 5, 3), 0.896825397, 1e-6);
    EXPECT_NEAR(stoch.leftCumulativeHypergeometric(10, 5, 5, 4), 0.996031746, 1e-6);
    EXPECT_NEAR(stoch.leftCumulativeHypergeometric(10, 5, 5, 5), 1, 1e-6);
}

TEST_F(CleanupGTest, testStochasticDist) {
    StochasticDistributionCalculator stoch(100);
    count kTotal = 10;
    count kIn = 3;
    count cOut = 20;
    count extStubs = 100;
    count M = extStubs - kTotal;
    auto fct = [](count x) {
        double f = 1.0;
        for (count i = 2; i <= x; ++i)
            f *= i;
        return f;
    };
    auto p = [&](count kIn) {
        count kOut = kTotal - kIn;
        count MInEdges = (M - cOut - kOut + kIn) / 2;
        double dividend = std::pow(2, -(double) kIn);
        double divisor = fct(kOut) * fct(kIn) * fct(cOut - kIn) * fct(MInEdges);
        return dividend / divisor;
    };
    double probabilitySum = 0.0;
    double cumulativeProbSum = 0.0;
    for (count x = 0; x < kTotal; ++x) {
        double probability = p(x);
        probabilitySum += probability;
        if (x >= kIn)
            cumulativeProbSum += probability;
    }
    double exactProbCorrect = p(kIn) / probabilitySum;
    double cumulativeProbCorrect = cumulativeProbSum / probabilitySum;

    double exactProb, cumulativeProb;
    std::tie(exactProb, cumulativeProb) = stoch.rightCumulativeStochastic(
            kTotal, kIn, cOut, extStubs);

    EXPECT_NEAR(exactProb, exactProbCorrect, 1e-6);
    EXPECT_NEAR(exactProb / exactProbCorrect, 1.0, 1e-6);
    EXPECT_NEAR(cumulativeProb, cumulativeProbCorrect, 1e-6);
    EXPECT_NEAR(cumulativeProb / cumulativeProbCorrect, 1.0, 1e-6);
}

} /* namespace NetworKit */
