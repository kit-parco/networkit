/*
 * CurveballGlobalTradeGenerator.cpp
 *
 *  Created on: Jul 11, 2017
 *      Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */

#include <algorithm>
#include <cassert>
#include <vector>

#include <networkit/randomization/CurveballGlobalTradeGenerator.hpp>

#include <networkit/auxiliary/Random.hpp>

namespace NetworKit {

CurveballGlobalTradeGenerator::value_type CurveballGlobalTradeGenerator::generate() const {
    using node_vector = std::vector<node>;

    value_type trades_out;
    trades_out.reserve(numGlobalTrades * numNodes / 2);

    node_vector node_permutation(numNodes);
    for (node node_id = 0; node_id < numNodes; node_id++) {
        node_permutation[node_id] = node_id;
    }

    for (count run = 0; run < numGlobalTrades; run++) {
        // shuffling a shuffled node_permutation is okay, no need
        // to reinitialize it
        std::shuffle(node_permutation.begin(), node_permutation.end(), Aux::Random::getURNG());

        auto rand_node_iter = node_permutation.cbegin();
        for (count t_id = 0; t_id < numNodes / 2; t_id++) {
            assert(rand_node_iter != node_permutation.cend());

            const node fst = *rand_node_iter;
            rand_node_iter++;
            const node snd = *rand_node_iter;
            rand_node_iter++;

            trades_out.emplace_back(fst, snd);
        }
    }

    return trades_out;
}

} // namespace NetworKit
