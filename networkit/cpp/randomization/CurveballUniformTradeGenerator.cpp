/*
 * CurveballUniformTradeGenerator.cpp
 *
 *  Created on: Jul 11, 2017
 *    Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */

#include <algorithm>
#include <cassert>
#include <vector>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/randomization/CurveballUniformTradeGenerator.hpp>

namespace NetworKit {

CurveballUniformTradeGenerator::value_type CurveballUniformTradeGenerator::generate() const {
    value_type trades_out;
    trades_out.clear();
    trades_out.reserve(numTrades);

    for (count t_id = 0; t_id < numTrades; t_id++) {
        const node fst = Aux::Random::integer(numNodes - 1);
        while (true) {
            const node snd = Aux::Random::integer(numNodes - 1);
            if (fst != snd) {
                trades_out.emplace_back(fst, snd);
                break;
            }
        }
    }

    return trades_out;
}

} // namespace NetworKit
