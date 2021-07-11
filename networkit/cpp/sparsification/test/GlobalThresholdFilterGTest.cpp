
#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/sparsification/GlobalThresholdFilter.hpp>
#include <networkit/sparsification/LocalDegreeScore.hpp>

namespace NetworKit {

class GlobalThresholdFilterGTest : public testing::TestWithParam<std::pair<bool, bool>> {
protected:
    bool isWeighted() const noexcept { return GetParam().first; }
    bool isDirected() const noexcept { return GetParam().second; }
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, GlobalThresholdFilterGTest,
                         testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                                         std::make_pair(false, true), std::make_pair(true, true)));

TEST_P(GlobalThresholdFilterGTest, testGlobalThresholdFilter) {
    static constexpr count n = 100;
    static constexpr double p = 0.15;

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(n, p, isDirected()).generate();
        if (isWeighted()) {
            G = GraphTools::toWeighted(G);
            G.forEdges([&G](const node u, const node v) {
                G.setWeight(u, v, Aux::Random::probability());
            });
        }

        G.indexEdges();

        LocalDegreeScore lds(G);
        lds.run();
        const auto scores = lds.scores();

        GlobalThresholdFilter gtf(G, scores, 0.5, false);
        auto G1 = gtf.calculate();

        EXPECT_EQ(G.isDirected(), G1.isDirected());
        EXPECT_EQ(G.isWeighted(), G1.isWeighted());
        EXPECT_EQ(G.numberOfNodes(), G1.numberOfNodes());
        EXPECT_EQ(G.upperNodeIdBound(), G1.upperNodeIdBound());

        G1.forEdges([&G](const node u, const node v, const edgeweight ew, edgeid) {
            EXPECT_TRUE(G.hasEdge(u, v));
            EXPECT_DOUBLE_EQ(G.weight(u, v), ew);
        });
    }
}

} // namespace NetworKit
