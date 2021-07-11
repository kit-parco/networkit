/*
 * CurveballUniformTradeGeneratorGTest.h
 *
 *  Created on: 26.05.2018
 *      Author:  Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */

#include <gtest/gtest.h>

#include <networkit/randomization/CurveballUniformTradeGenerator.hpp>

namespace NetworKit {

class CurveballUniformTradeGeneratorGTest : public ::testing::Test {};

TEST_F(CurveballUniformTradeGeneratorGTest, testGeneration) {
    CurveballUniformTradeGenerator gen(15, 10);
    auto trades = gen.generate();
    ASSERT_EQ(trades.size(), 15u);
    for (auto t_it = trades.cbegin(); t_it != trades.cend(); t_it++) {
        ASSERT_LE((*t_it).first, 9u);
        ASSERT_GE((*t_it).first, 0u);
        ASSERT_LE((*t_it).second, 9u);
        ASSERT_GE((*t_it).second, 0u);
    }
}

} // namespace NetworKit
