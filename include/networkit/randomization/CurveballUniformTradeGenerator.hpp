/*
 * CurveballUniformTradeGenerator.hpp
 *
 *  Created on: Jul 11, 2017
 *      Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */

#ifndef NETWORKIT_RANDOMIZATION_CURVEBALL_UNIFORM_TRADE_GENERATOR_HPP_
#define NETWORKIT_RANDOMIZATION_CURVEBALL_UNIFORM_TRADE_GENERATOR_HPP_

#include <stdexcept>
#include <utility>
#include <vector>

#include <networkit/Globals.hpp>

namespace NetworKit {

class CurveballUniformTradeGenerator final {
public:
    using value_type = std::vector<std::pair<node, node>>;

private:
    const count numTrades;
    const node numNodes;

public:
    CurveballUniformTradeGenerator(count trade_num, node num_nodes)
        : numTrades(trade_num), numNodes(num_nodes) {
        if (num_nodes <= 1)
            throw std::runtime_error("At least two nodes are required.");
    }

    value_type generate() const;
};

} // namespace NetworKit

#endif // NETWORKIT_RANDOMIZATION_CURVEBALL_UNIFORM_TRADE_GENERATOR_HPP_
