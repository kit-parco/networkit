// no-networkit-format
/*
 * MatcherGTest.cpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Henning
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/matching/Matcher.hpp>
#include <networkit/matching/Matching.hpp>
#include <networkit/matching/PathGrowingMatcher.hpp>
#include <networkit/matching/LocalMaxMatcher.hpp>
#include <networkit/matching/SuitorMatcher.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/DibapGraphReader.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class MatcherGTest: public testing::Test {};

TEST_F(MatcherGTest, testLocalMaxMatching) {
    {
        Graph G(10, true, true);
        EXPECT_THROW(LocalMaxMatcher{G}, std::runtime_error);
    }

    count n = 50;
    Graph G(n);
    G.forNodePairs([&](node u, node v){
        G.addEdge(u,v);
    });

    LocalMaxMatcher localMaxMatcher(G);

    TRACE("Start localMax matching");
    localMaxMatcher.run();
    Matching M = localMaxMatcher.getMatching();
    TRACE("Finished localMax matching");

    count numExpEdges = n / 2;
    bool isProper = M.isProper(G);
    EXPECT_TRUE(isProper);
    EXPECT_EQ(M.size(G), numExpEdges);

#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
    DibapGraphReader reader;
    Graph airfoil1 = reader.read("input/airfoil1.gi");
    LocalMaxMatcher lmm(airfoil1);
    lmm.run();
    M = lmm.getMatching();
    isProper = M.isProper(airfoil1);
    EXPECT_TRUE(isProper);
    DEBUG("LocalMax on airfoil1 produces matching of size: " , M.size(G));
#endif
}

TEST_F(MatcherGTest, testLocalMaxMatchingDirectedWarning) {
    Graph G(2, false, true);
    G.addEdge(0,1);
    EXPECT_THROW(LocalMaxMatcher localMaxMatcher(G), std::runtime_error);
}

TEST_F(MatcherGTest, testPgaMatchingOnWeightedGraph) {
    count n = 50;
    Graph G(n);
    G.forNodePairs([&](node u, node v){
        G.addEdge(u,v, Aux::Random::real());
    });
    PathGrowingMatcher pgaMatcher(G);
    EXPECT_NO_THROW(pgaMatcher.run());
}

TEST_F(MatcherGTest, testPgaMatchingWithSelfLoops) {
    count n = 50;
    Graph G(n);
    G.forNodePairs([&](node u, node v){
        G.addEdge(u,v, Aux::Random::real());
    });
    G.forNodes([&](node u){
        G.addEdge(u,u);
    });
    EXPECT_THROW(PathGrowingMatcher pgaMatcher(G),std::invalid_argument);
}


TEST_F(MatcherGTest, testPgaMatching) {
    count n = 50;
    Graph G(n);
    G.forNodePairs([&](node u, node v){
        G.addEdge(u,v);
    });
    PathGrowingMatcher pgaMatcher(G);

    DEBUG("Start PGA matching on 50-clique");

    pgaMatcher.run();
    Matching M = pgaMatcher.getMatching();

    count numExpEdges = n / 2;
    bool isProper = M.isProper(G);
    EXPECT_TRUE(isProper);
    EXPECT_EQ(M.size(G), numExpEdges);
    DEBUG("Finished PGA matching on 50-clique");


#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
    DibapGraphReader reader;
    Graph airfoil1 = reader.read("input/airfoil1.gi");
    PathGrowingMatcher pga2(airfoil1);
    pga2.run();
    M = pga2.getMatching();
    isProper = M.isProper(airfoil1);
    EXPECT_TRUE(isProper);
    DEBUG("PGA on airfoil1 produces matching of size: " , M.size(G));
#endif
}

TEST_F(MatcherGTest, debugValidMatching) {
    METISGraphReader reader;
    Graph G = reader.read("coAuthorsDBLP.graph");

    LocalMaxMatcher pmatcher(G);
    pmatcher.run();
    Matching M = pmatcher.getMatching();

    bool isProper = M.isProper(G);
    EXPECT_TRUE(isProper);
}

TEST_F(MatcherGTest, testSuitorMatcher) {
    {   // Directed graphs are not supported
        Graph G(10, true, true);
        EXPECT_THROW(SuitorMatcher{G}, std::runtime_error);
    }

    {   // Graphs with self loops are not supported
        Graph G(10);
        G.addEdge(0, 0);
        G.addEdge(0, 0);
        EXPECT_THROW(SuitorMatcher{G}, std::runtime_error);
    }

    const edgeweight maxWeight = 10;

    const auto hasUnmatchedNeighbors = [](const Graph &G, const Matching &M) -> bool {
        for (const auto e : G.edgeRange())
            if (!M.isMatched(e.u) && !M.isMatched(e.v))
                return true;
        return false;
    };

    const auto doTest = [maxWeight, &hasUnmatchedNeighbors](Graph &G) -> void {
        // Test suitor matcher
        SuitorMatcher sm(G, false, true);
        sm.run();
        const auto M1 = sm.getMatching();
        EXPECT_TRUE(M1.isProper(G));
        EXPECT_FALSE(hasUnmatchedNeighbors(G, M1));

        GraphTools::sortEdgesByWeight(G, true);
        {
            auto G1 = G;
            G1.addEdge(0, 1, maxWeight);
            G1.addEdge(0, 2, maxWeight);
            EXPECT_THROW(SuitorMatcher(G1, true, true), std::runtime_error);
        }

        // Test sort suitor matcher
        SuitorMatcher ssm(G, true, true);
        ssm.run();
        const auto M2 = ssm.getMatching();
        EXPECT_TRUE(M2.isProper(G));
        EXPECT_FALSE(hasUnmatchedNeighbors(G, M2));

        // Matchings must be the same
        G.forNodes([&M1, &M2](node u) { EXPECT_EQ(M1.mate(u), M2.mate(u)); });
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);

        // Test unweighted
        auto G = METISGraphReader{}.read("input/PGPgiantcompo.graph");
        G.removeSelfLoops();
        G.removeMultiEdges();
        doTest(G);

        // Test weighted
        G = GraphTools::toWeighted(G);
        G.forEdges([&G, maxWeight](node u, node v) {
            G.setWeight(u, v, Aux::Random::real(maxWeight));
        });
        doTest(G);
    }
}

} // namespace NetworKit
