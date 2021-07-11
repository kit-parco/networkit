/*
 * CurveballImplGTest.cpp
 *
 *  Created on: 26.08.2018
 *      Author:  Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */

#include <gtest/gtest.h>

#include "../CurveballImpl.hpp"

namespace NetworKit {

class CurveballImplGTest : public ::testing::Test {};

namespace CurveballDetails {

using degree_vector = std::vector<count>;
using neighbour_it = std::vector<node>::const_iterator;

TEST_F(CurveballImplGTest, testContainer) {
    degree_vector degrees = {3, 2, 1, 1, 1};

    const edgeid degree_sum = 8;

    // =========== Container ===========

    CurveballAdjacencyList adj_list(degrees, degree_sum);

    /*
     * Expected Structure:
     * | 0           | 1        | 2     | 3     | 4     |
     * +-------------+----------+-------+-------+-------+
     * | _, _, _, END| _, _, END| _, END| _, END| _, END|
     */

    // Node 0
    neighbour_it n_it = adj_list.cbegin(0);
    ++n_it;
    ++n_it;
    ++n_it;
    ASSERT_EQ(*n_it, LISTROW_END);
    n_it++;

    // Node 1
    n_it++;
    n_it++;
    ASSERT_EQ(*n_it, LISTROW_END);
    n_it++;

    // Node 2
    n_it++;
    ASSERT_EQ(*n_it, LISTROW_END);
    n_it++;

    // Node 3
    n_it++;
    ASSERT_EQ(*n_it, LISTROW_END);
    n_it++;

    // Node 4
    n_it++;
    ASSERT_EQ(*n_it, LISTROW_END);

    // =========== Insertion ===========

    adj_list.insertNeighbour(0, 3);
    adj_list.insertNeighbour(0, 2);

    /*
     * Expected Structure:
     * | 0           | 1        | 2     | 3     | 4     |
     * +-------------+----------+-------+-------+-------+
     * | 3, 2, _, END| _, _, END| _, END| _, END| _, END|
     */

    n_it = adj_list.cbegin(0);
    ASSERT_EQ(*n_it, 3);
    n_it++;
    ASSERT_EQ(*n_it, 2);
    n_it++;
    ASSERT_EQ(n_it, adj_list.cend(0));

    // =========== Reset ===========

    adj_list.resetRow(0);
    ASSERT_EQ(adj_list.cbegin(0), adj_list.cend(0));
}

TEST_F(CurveballImplGTest, testContainerByInitialize) {
    degree_vector degrees = {3, 2, 1, 1, 1};

    const edgeid degree_sum = 8;

    // =========== Container ===========

    CurveballAdjacencyList adj_list;
    adj_list.initialize(degrees, degree_sum);

    /*
     * Expected Structure:
     * | 0           | 1        | 2     | 3     | 4     |
     * +-------------+----------+-------+-------+-------+
     * | _, _, _, END| _, _, END| _, END| _, END| _, END|
     */

    // Node 0
    neighbour_it n_it = adj_list.cbegin(0);
    ++n_it;
    ++n_it;
    ++n_it;
    ASSERT_EQ(*n_it, LISTROW_END);
    n_it++;

    // Node 1
    n_it++;
    n_it++;
    ASSERT_EQ(*n_it, LISTROW_END);
    n_it++;

    // Node 2
    n_it++;
    ASSERT_EQ(*n_it, LISTROW_END);
    n_it++;

    // Node 3
    n_it++;
    ASSERT_EQ(*n_it, LISTROW_END);
    n_it++;

    // Node 4
    n_it++;
    ASSERT_EQ(*n_it, LISTROW_END);

    // =========== Insertion ===========

    adj_list.insertNeighbour(0, 3);
    adj_list.insertNeighbour(0, 2);

    /*
     * Expected Structure:
     * | 0           | 1        | 2     | 3     | 4     |
     * +-------------+----------+-------+-------+-------+
     * | 3, 2, _, END| _, _, END| _, END| _, END| _, END|
     */

    n_it = adj_list.cbegin(0);
    ASSERT_EQ(*n_it, 3);
    n_it++;
    ASSERT_EQ(*n_it, 2);
    n_it++;
    ASSERT_EQ(n_it, adj_list.cend(0));

    // =========== Reset ===========

    adj_list.resetRow(0);
    ASSERT_EQ(adj_list.cbegin(0), adj_list.cend(0));
}

///////////////////////////////////////////////////////////////////////////////

TEST_F(CurveballImplGTest, testTradeListConstructor) {
    std::vector<trade_descriptor> trades = {{1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 5}, {5, 8}};

    /*
     * Expected Structure:
     * Node     Trade-ids it partakes in
     * 0:       END
     * 1:       0, 1, 2, 3, END
     * 2:       0, 4, END
     * 3:       1, END
     * 4:       2, END
     * 5:       3, 4, 5, END
     * 6:       END
     * 7:       END
     * 8:       5, END
     * 9:       END
     */

    // more nodes
    const node num_nodes = 10;

    TradeList trade_list(trades, num_nodes);

    auto trade_iter = trade_list.getTrades(0);
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    trade_iter = trade_list.getTrades(1);
    ASSERT_EQ(*trade_iter, 0);
    trade_iter++;
    ASSERT_EQ(*trade_iter, 1);
    trade_iter++;
    ASSERT_EQ(*trade_iter, 2);
    trade_iter++;
    ASSERT_EQ(*trade_iter, 3);
    trade_iter++;
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    trade_iter = trade_list.getTrades(2);
    ASSERT_EQ(*trade_iter, 0);
    trade_iter++;
    ASSERT_EQ(*trade_iter, 4);
    trade_iter++;
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    trade_iter = trade_list.getTrades(3);
    ASSERT_EQ(*trade_iter, 1);
    trade_iter++;
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    trade_iter = trade_list.getTrades(4);
    ASSERT_EQ(*trade_iter, 2);
    trade_iter++;
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    trade_iter = trade_list.getTrades(5);
    ASSERT_EQ(*trade_iter, 3);
    trade_iter++;
    ASSERT_EQ(*trade_iter, 4);
    trade_iter++;
    ASSERT_EQ(*trade_iter, 5);
    trade_iter++;
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    ASSERT_EQ(*(trade_list.getTrades(6)), TRADELIST_END);

    ASSERT_EQ(*(trade_list.getTrades(7)), TRADELIST_END);

    trade_iter = trade_list.getTrades(8);
    ASSERT_EQ(*trade_iter, 5);
    trade_iter++;
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    ASSERT_EQ(*(trade_list.getTrades(9)), TRADELIST_END);
}

TEST_F(CurveballImplGTest, testTradeListInitialize) {
    std::vector<trade_descriptor> trades = {{1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 5}, {5, 8}};

    /*
     * Expected Structure:
     * Node     Trade-ids it partakes in
     * 0:       END
     * 1:       0, 1, 2, 3, END
     * 2:       0, 4, END
     * 3:       1, END
     * 4:       2, END
     * 5:       3, 4, 5, END
     * 6:       END
     * 7:       END
     * 8:       5, END
     * 9:       END
     */

    // more nodes
    const node num_nodes = 10;

    TradeList trade_list(num_nodes);
    trade_list.initialize(trades);

    auto trade_iter = trade_list.getTrades(0);
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    trade_iter = trade_list.getTrades(1);
    ASSERT_EQ(*trade_iter, 0);
    trade_iter++;
    ASSERT_EQ(*trade_iter, 1);
    trade_iter++;
    ASSERT_EQ(*trade_iter, 2);
    trade_iter++;
    ASSERT_EQ(*trade_iter, 3);
    trade_iter++;
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    trade_iter = trade_list.getTrades(2);
    ASSERT_EQ(*trade_iter, 0);
    trade_iter++;
    ASSERT_EQ(*trade_iter, 4);
    trade_iter++;
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    trade_iter = trade_list.getTrades(3);
    ASSERT_EQ(*trade_iter, 1);
    trade_iter++;
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    trade_iter = trade_list.getTrades(4);
    ASSERT_EQ(*trade_iter, 2);
    trade_iter++;
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    trade_iter = trade_list.getTrades(5);
    ASSERT_EQ(*trade_iter, 3);
    trade_iter++;
    ASSERT_EQ(*trade_iter, 4);
    trade_iter++;
    ASSERT_EQ(*trade_iter, 5);
    trade_iter++;
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    ASSERT_EQ(*(trade_list.getTrades(6)), TRADELIST_END);

    ASSERT_EQ(*(trade_list.getTrades(7)), TRADELIST_END);

    trade_iter = trade_list.getTrades(8);
    ASSERT_EQ(*trade_iter, 5);
    trade_iter++;
    ASSERT_EQ(*trade_iter, TRADELIST_END);

    ASSERT_EQ(*(trade_list.getTrades(9)), TRADELIST_END);
}

} // namespace CurveballDetails
} // namespace NetworKit
