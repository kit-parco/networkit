
#include <gtest/gtest.h>

#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/reachability/ReachableNodes.hpp>

namespace NetworKit {

class ReachabilityGTest : public testing::TestWithParam<std::pair<bool, bool>> {
protected:
    bool isDirected() const noexcept;
    bool isWeighted() const noexcept;
    static void generateRandomWeights(Graph &G);
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, ReachabilityGTest,
                         testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                                         std::make_pair(false, true), std::make_pair(true, true)));

bool ReachabilityGTest::isWeighted() const noexcept {
    return GetParam().first;
}

bool ReachabilityGTest::isDirected() const noexcept {
    return GetParam().second;
}

void ReachabilityGTest::generateRandomWeights(Graph &G) {
    if (!G.isWeighted())
        G = GraphTools::toWeighted(G);
    G.forEdges([&G = G](node u, node v) { G.setWeight(u, v, Aux::Random::real(1, 10)); });
}

TEST_P(ReachabilityGTest, testReachableNodesExact) {
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(100, 0.05, isDirected()).generate();
        if (isWeighted())
            generateRandomWeights(G);

        ReachableNodes rb(G);
        rb.run();

        G.balancedParallelForNodes([&](node u) {
            count r = 0;
            Traversal::BFSfrom(G, u, [&r = r](node, count) { ++r; });
            EXPECT_EQ(rb.numberOfReachableNodes(u), r);
        });
    }
}

TEST_P(ReachabilityGTest, testReachableNodesApprox) {
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(100, 0.05, isDirected()).generate();
        if (isWeighted())
            generateRandomWeights(G);

        ReachableNodes rb(G, false);
        rb.run();

        G.balancedParallelForNodes([&](node u) {
            count r = 0;
            Traversal::BFSfrom(G, u, [&r = r](node, count) { ++r; });
            EXPECT_LE(rb.numberOfReachableNodesLB(u), r);
            EXPECT_GE(rb.numberOfReachableNodesUB(u), r);
        });
    }
}

} // namespace NetworKit
