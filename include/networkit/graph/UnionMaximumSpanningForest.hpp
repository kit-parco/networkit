// no-networkit-format
#ifndef NETWORKIT_GRAPH_UNION_MAXIMUM_SPANNING_FOREST_HPP_
#define NETWORKIT_GRAPH_UNION_MAXIMUM_SPANNING_FOREST_HPP_

#include <limits>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/UnionFind.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * Union maximum-weight spanning forest algorithm, computes the union of all maximum-weight spanning forests using Kruskal's algorithm.
 */
class UnionMaximumSpanningForest final : public Algorithm {
public:
    /**
     * Initialize the union maximum-weight spanning forest algorithm, uses edge weights.
     *
     * @param G The input graph.
     */
    UnionMaximumSpanningForest(const Graph &G);

    /**
     * Initialize the union maximum-weight spanning forest algorithm using an attribute as edge weight.
     *
     * This copies the attribute values, the supplied attribute vector is not stored.
     *
     * @param G The input graph.
     * @param attribute The attribute to use, can be either of type edgeweight (double) or count (uint64), internally all values are handled as double.
     */
    template <typename A>
    UnionMaximumSpanningForest(const Graph &G, const std::vector<A> &attribute);

    /**
     * Execute the algorithm. The algorithm is not parallel.
     */
    void run() override;

    /**
     * Get a boolean attribute that indicates for each edge if it is part of any maximum-weight spanning forest.
     *
     * This attribute is only calculated and can thus only be request if the supplied graph has edge ids.
     *
     * @param move If the attribute shall be moved out of the algorithm instance.
     * @return The vector with the boolean attribute for each edge.
     */
    std::vector<bool> getAttribute(bool move = false);

    /**
     * Checks if the edge (@a u, @a v) is part of any maximum-weight spanning forest.
     *
     * @param u The first node of the edge to check
     * @param v The second node of the edge to check
     * @return If the edge is part of any maximum-weight spanning forest.
     */
    bool inUMSF(node u, node v) const;

    /**
     * Checks if the edge with the id @a eid is part of any maximum-weight spanning forest.
     *
     * @param eid The id of the edge to check.
     * @return If the edge is part of any maximum-weight spanning forest.
     */
    bool inUMSF(edgeid eid) const;

    /**
     * Gets the union of all maximum-weight spanning forests as graph.
     *
     * @param move If the graph shall be moved out of the algorithm instance.
     * @return The calculated union of all maximum-weight spanning forests.
     */
    Graph getUMSF(bool move = false);

    /**
     * @return false - this algorithm is not parallelized.
     */
    bool TLX_DEPRECATED(isParallel() const override);

    /**
     * @return The name of the algorithm.
     */
    std::string TLX_DEPRECATED(toString() const override);

private:
    struct weightedEdge {
        edgeweight attribute;
        node u;
        node v;
        edgeid eid;

        bool operator>(const weightedEdge &other) const {
            return (attribute > other.attribute) || (attribute == other.attribute && (u > other.u || (u == other.u && v > other.v)));
        };
        weightedEdge(node u, node v, edgeweight attribute, edgeid eid = 0) : attribute(attribute), u(u), v(v), eid(eid) {};
    };

    const Graph *G;
    std::vector<weightedEdge> weightedEdges;

    Graph umsf;
    std::vector<bool> umsfAttribute;

    bool hasWeightedEdges;
    bool hasUMSF;
    bool hasAttribute;
};

template <typename A>
UnionMaximumSpanningForest::UnionMaximumSpanningForest(const Graph &G, const std::vector< A > &attribute) : G(&G), hasWeightedEdges(false), hasUMSF(false), hasAttribute(false) {
    if (!G.hasEdgeIds()) {
        throw std::runtime_error("Error: Edges of G must be indexed for using edge attributes");
    }

    weightedEdges.reserve(G.numberOfEdges());

    G.forEdges([&](node u, node v, edgeid eid) {
        weightedEdges.emplace_back(u, v, attribute[eid], eid);
    });

    INFO(weightedEdges.size(), " weighted edges saved");

    hasWeightedEdges = true;
}

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_UNION_MAXIMUM_SPANNING_FOREST_HPP_
