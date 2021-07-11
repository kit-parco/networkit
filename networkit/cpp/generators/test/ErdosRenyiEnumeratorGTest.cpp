// no-networkit-format
/*
 * ErdosReyniEnumeratorGTest.cpp
 *
 *  Created on: 09. Aug. 2018
 *      Author: Manuel Penschuck
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <numeric>
#include <tuple>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ErdosRenyiEnumerator.hpp>

namespace NetworKit {

class ErdosRenyiEnumeratorGTest : public testing::TestWithParam<std::tuple<bool, node, double> > {};

template <bool Parallel, bool FixedPoint>
static void testEre(const bool directed, const node n, const double prob) {
    Aux::Random::setSeed(static_cast<size_t>((123456 + directed * 23425 + 234*Parallel + 34*FixedPoint) * prob * n), false);

    std::vector<unsigned> edge_counts(n * n, 0);
    std::vector<unsigned> edge_exists(n * n, 0); // std::vector<bool> is not thread safe

    auto edgeindex = [](node u, node v, node n) -> size_t {
        return u*n + v;
    };

    // after this many rounds each edge should be hit at least once
    const auto rounds = (0.0 < prob && prob < 1.0) ? static_cast<int>(20.0 / prob) : 1;

    for(int i=0; i<rounds; i++) {
        ErdosRenyiEnumerator<FixedPoint> ere(n, prob, directed);

        // initialize data
        edge_exists.assign(n*n, 0);

        constexpr size_t cacheline_scale = 17; // avoid false sharing
        std::vector<count> num_edges_thread(omp_get_max_threads() * cacheline_scale);

        auto handle = [&] (int tid, node u, node v) {
            // ensure edge is legal (for undirected edges the
            // first node is larger than the second)
            ASSERT_LE(u, n);
            ASSERT_LE(v, directed ? n : u);

            // check that edge is unique
            ASSERT_FALSE(edge_exists[edgeindex(u,v,n)]);
            edge_exists[edgeindex(u,v,n)] = 1;

            // collect statistics for probabilisitic checks
            edge_counts[edgeindex(u,v,n)]++;
            num_edges_thread[cacheline_scale*tid]++;
        };

        count num_edges_gen;
        if (Parallel)
            num_edges_gen = ere.forEdgesParallel(handle);
        else
            num_edges_gen = ere.forEdges(handle);

        size_t num_edges = std::accumulate(num_edges_thread.cbegin(), num_edges_thread.cend(), 0u);
        size_t active_threads = std::count_if(num_edges_thread.cbegin(), num_edges_thread.cend(), [] (count c) {return !!c;});

        ASSERT_EQ(num_edges, num_edges_gen);

        // Check that result is somewhat balanced along threads
        if (Parallel && prob > 0) {
            EXPECT_EQ(omp_get_max_threads(), static_cast<int>(active_threads));
            for (auto count : num_edges_thread)
                ASSERT_LE(count, 2 * num_edges / active_threads);
        }

        if (prob == 0.0 || prob == 1.0) {
            ASSERT_EQ(num_edges, ere.expectedNumberOfEdges());
        } else {
            // count that we produced a number of edges, near the expected number
            // we could apply much sharper bounds here if we allow for rare errors
            ASSERT_NEAR(num_edges, ere.expectedNumberOfEdges(), 0.6 * ere.expectedNumberOfEdges());
        }
    }

    if (prob == 0.0) return;

    for(node u = 1; u != n; u++) {
        const node upper = directed ? n : u;
        for(node v = 0; v != upper; v++) {
            // Again, much sharper bounds are possible if we do not
            // rely on the implicit union bound here.
            auto count = edge_counts[edgeindex(u,v,n)];
            EXPECT_GE(count, 1u) << "edge(" << u << ", " << v << "), rounds=" << rounds;
            if (prob < 1.0) {
                ASSERT_LE(count, log(rounds) * rounds * prob) << "edge(" << u << ", " << v << "), rounds=" << rounds;
            }
        }
    }
}

TEST_P(ErdosRenyiEnumeratorGTest, TestFloatingSequential) {
    testEre<false, false>(std::get<0>(GetParam()), std::get<1>(GetParam()), std::get<2>(GetParam()));
}

TEST_P(ErdosRenyiEnumeratorGTest, TestFloatingParallel) {
    testEre<true, false>(std::get<0>(GetParam()), std::get<1>(GetParam()), std::get<2>(GetParam()));
}

TEST_P(ErdosRenyiEnumeratorGTest, TestFixedPointSequential) {
    testEre<false, true>(std::get<0>(GetParam()), std::get<1>(GetParam()), std::get<2>(GetParam()));
}

TEST_P(ErdosRenyiEnumeratorGTest, TestFixedPointParallel) {
    testEre<true, true>(std::get<0>(GetParam()), std::get<1>(GetParam()), std::get<2>(GetParam()));
}


INSTANTIATE_TEST_SUITE_P(ErdosRenyiEnumeratorGTest, ErdosRenyiEnumeratorGTest,
                        ::testing::Values(
                         std::make_tuple(false, 100, 0.0),  std::make_tuple(true, 100, 0.0),
                         std::make_tuple(false, 100, 0.1),  std::make_tuple(true, 100, 0.1),
                         std::make_tuple(false, 100, 0.5),  std::make_tuple(true, 100, 0.5),
                         std::make_tuple(false, 100, 0.7),  std::make_tuple(true, 100, 0.7),
                         std::make_tuple(false, 100, 1.0),  std::make_tuple(true, 100, 1.0),
                         std::make_tuple(false, 200, 0.01), std::make_tuple(true, 200, 0.01)
                        ));

} // ! namespace NetworKit
