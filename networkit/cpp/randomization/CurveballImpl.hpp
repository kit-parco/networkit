/*
 * CurveballImpl.hpp
 *
 * Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */
#ifndef RANDOMIZATION_CURVEBALL_IMPL_H
#define RANDOMIZATION_CURVEBALL_IMPL_H

#include <cassert>
#include <utility>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace CurveballDetails {

// Global Definitions
using trade_descriptor = std::pair<node, node>;
using tradeid = node;
using count = node;

using trade_vector = std::vector<trade_descriptor>;
using nodepair_vector = std::vector<std::pair<node, node>>;

constexpr node INVALID_NODE = std::numeric_limits<node>::max();
constexpr count LISTROW_END = std::numeric_limits<count>::max();
constexpr tradeid TRADELIST_END = std::numeric_limits<tradeid>::max();

class CurveballAdjacencyList {
public:
    using degree_vector = std::vector<count>;
    using neighbour_vector = std::vector<node>;
    using pos_vector = std::vector<edgeid>;
    using pos_it = pos_vector::iterator;
    using neighbour_it = neighbour_vector::iterator;
    using cneighbour_it = neighbour_vector::const_iterator;
    using nodepair_vector = std::vector<std::pair<node, node>>;

private:
    neighbour_vector neighbours;
    degree_vector offsets;
    pos_vector begins;
    edgeid degreeCount;

public:
    CurveballAdjacencyList() = default;

    // Receives the degree_vector to initialize
    // As trades permute neighbours the degrees don't change
    CurveballAdjacencyList(const degree_vector &degrees, edgeid degree_count);

    void initialize(const degree_vector &degrees, edgeid degree_count);

    void restructure();

    // No Copy Constructor
    CurveballAdjacencyList(const CurveballAdjacencyList &) = delete;

    neighbour_it begin(node node_id) { return neighbours.begin() + begins[node_id]; }

    neighbour_it end(node node_id) {
        return neighbours.begin() + begins[node_id] + offsets[node_id];
    }

    cneighbour_it cbegin(node node_id) const { return neighbours.cbegin() + begins[node_id]; }

    cneighbour_it cend(node node_id) const {
        return neighbours.cbegin() + begins[node_id] + offsets[node_id];
    }

    nodepair_vector getEdges() const;

    void insertNeighbour(node node_id, node neighbour) {
        auto pos = begin(node_id) + offsets[node_id];

        assert(*pos != LISTROW_END);

        *pos = neighbour;

        offsets[node_id]++;
    }

    node numberOfNodes() const { return static_cast<node>(offsets.size()); }

    node numberOfEdges() const { return static_cast<edgeid>(degreeCount); }

    void resetRow(node node_id) {
        assert(node_id < static_cast<node>(offsets.size()));

        offsets[node_id] = 0;

        return;
    }

    count degreeAt(node node_id) const {
        assert(node_id < static_cast<node>(offsets.size()));

        return begins[node_id + 1] - begins[node_id] - 1;
    }
};

class CurveballMaterialization {

private:
    const CurveballAdjacencyList &adjacencyList;

public:
    CurveballMaterialization(const CurveballAdjacencyList &adj_list);

    Graph toGraph(bool parallel);

private:
    void toGraphParallel(Graph &G);
    void toGraphSequential(Graph &G);
};

class TradeList {
public:
    using edge_vector = std::vector<std::pair<node, node>>;
    using offset_vector = std::vector<tradeid>;
    using tradeid_vector = std::vector<tradeid>;
    using trade = trade_descriptor;
    using trade_vector = std::vector<trade>;
    using tradeid_it = tradeid_vector::const_iterator;

private:
    tradeid_vector tradeList;
    offset_vector offsets;
    const node numNodes;

public:
    TradeList(node num_nodes);

    // Receives the edge_vector to initialize
    TradeList(const trade_vector &trades, node num_nodes);

    // Initialize method
    void initialize(const trade_vector &trades);

    // No Copy Constructor
    TradeList(const TradeList &) = delete;

    tradeid_it getTrades(node nodeid) const {
        assert(nodeid < numNodes);

        return tradeList.begin() + offsets[nodeid];
    }

    void incrementOffset(node nodeid) {
        assert(nodeid < numNodes);
        assert(1 <= offsets[nodeid + 1] - offsets[nodeid]);

        offsets[nodeid]++;
    }

    node numberOfNodes() const { return numNodes; }
};

class CurveballIM {
public:
    CurveballIM(const Graph &G);

    void run(const trade_vector &trades);

    count getNumberOfAffectedEdges() const {
        assert(hasRun);
        return numAffectedEdges;
    }

    Graph getGraph(bool parallel) const;

    nodepair_vector getEdges() const;

private:
    const Graph *G;
    const node numNodes;

    bool hasRun;
    CurveballAdjacencyList adjList;
    TradeList tradeList;
    count maxDegree;
    edgeid numAffectedEdges; // affected half-edges

    void loadFromGraph(const trade_vector &trades);

    void restructureGraph(const trade_vector &trades);

    inline void update(node a, node b) {
        const tradeid ta = *(tradeList.getTrades(a));
        const tradeid tb = *(tradeList.getTrades(b));
        if (ta < tb) {
            adjList.insertNeighbour(a, b);
            return;
        }

        if (ta > tb) {
            adjList.insertNeighbour(b, a);
            return;
        }
        { adjList.insertNeighbour(a, b); }
    }
};

} // namespace CurveballDetails
} // namespace NetworKit

#endif // ! RANDOMIZATION_CURVEBALL_IMPL_H
