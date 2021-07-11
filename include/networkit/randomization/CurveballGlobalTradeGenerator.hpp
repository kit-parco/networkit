/*
 * CurveballGlobalTradeGenerator.hpp
 *
 *  Created on: Jul 11, 2017
 *      Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */

#ifndef NETWORKIT_RANDOMIZATION_CURVEBALL_GLOBAL_TRADE_GENERATOR_HPP_
#define NETWORKIT_RANDOMIZATION_CURVEBALL_GLOBAL_TRADE_GENERATOR_HPP_

#include <stdexcept>
#include <utility>
#include <vector>

#include <networkit/Globals.hpp>

namespace NetworKit {

class CurveballGlobalTradeGenerator final {
public:
    using value_type = std::vector<std::pair<node, node>>;

private:
    const node numNodes;
    const count numGlobalTrades;

public:
    CurveballGlobalTradeGenerator(count run_length, node num_nodes)
        : numNodes(num_nodes), numGlobalTrades(run_length) {
        if (num_nodes <= 1)
            throw std::runtime_error("At least two nodes are required.");
    }

    value_type generate() const;
};

} // namespace NetworKit

#endif // NETWORKIT_RANDOMIZATION_CURVEBALL_GLOBAL_TRADE_GENERATOR_HPP_
