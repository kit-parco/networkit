#ifndef NETWORKIT_GRAPH_RANDOM_MAXIMUM_SPANNING_FOREST_HPP_
#define NETWORKIT_GRAPH_RANDOM_MAXIMUM_SPANNING_FOREST_HPP_

#include <limits>

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/UnionFind.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/base/Algorithm.hpp>

namespace NetworKit {

/**
 * Computes a random maximum-weight spanning forest using Kruskal's algorithm by randomizing the order of edges of the same weight.
 */
class RandomMaximumSpanningForest : public Algorithm {
public:
    /**
     * Initialize the random maximum-weight spanning forest algorithm, uses edge weights.
     *
     * @param G The input graph.
     */
    RandomMaximumSpanningForest(const Graph &G);

    /**
     * Initialize the random maximum-weight spanning forest algorithm using an attribute as edge weight.
     *
     * This copies the attribute values, the supplied attribute vector is not stored.
     *
     * @param G The input graph.
     * @param attribute The attribute to use, can be either of type edgeweight (double) or count (uint64), internally all values are handled as double.
     */
    template <typename A>
    RandomMaximumSpanningForest(const Graph &G, const std::vector<A> &attribute);

    /**
     * Execute the algorithm.
     */
    virtual void run() override;

    /**
     * Get a boolean attribute that indicates for each edge if it is part of the calculated maximum-weight spanning forest.
     *
     * This attribute is only calculated and can thus only be request if the supplied graph has edge ids.
     *
     * @param move If the attribute shall be moved out of the algorithm instance.
     * @return The vector with the boolean attribute for each edge.
     */
    std::vector<bool> getAttribute(bool move = false);

    /**
     * Checks if the edge (@a u, @a v) is part of the calculated maximum-weight spanning forest.
     *
     * @param u The first node of the edge to check
     * @param v The second node of the edge to check
     * @return If the edge is part of the calculated maximum-weight spanning forest.
     */
    bool inMSF(node u, node v) const;

    /**
     * Checks if the edge with the id @a eid is part of the calculated maximum-weight spanning forest.
     *
     * @param eid The id of the edge to check.
     * @return If the edge is part of the calculated maximum-weight spanning forest.
     */
    bool inMSF(edgeid eid) const;

    /**
     * Gets the calculated maximum-weight spanning forest as graph.
     *
     * @return The calculated maximum-weight spanning forest.
     */
    const Graph& getMSF() const;

    /**
     * @return false - this algorithm is not parallelized
     */
    virtual bool isParallel() const override;

    /**
     * @return The name of this algorithm.
     */
    virtual std::string toString() const override;

private:
    struct weightedEdge {
        double attribute;
        node u;
        node v;
        edgeid eid;
        index rand;

        bool operator>(const weightedEdge &other) const {
            return
            (attribute > other.attribute)
            ||
            (attribute == other.attribute &&
                (rand > other.rand ||
                (rand == other.rand && (u > other.u || (u == other.u && v > other.v)))));
        };
        weightedEdge(node u, node v, double attribute, edgeid eid = 0) : attribute(attribute), u(u), v(v), eid(eid), rand(Aux::Random::integer()) {};
    };

    const Graph &G;
    std::vector<weightedEdge> weightedEdges;

    Graph msf;
    std::vector<bool> msfAttribute;

    bool hasWeightedEdges;
    bool hasMSF;
    bool hasAttribute;
};

template <typename A>
RandomMaximumSpanningForest::RandomMaximumSpanningForest(const Graph &G, const std::vector< A > &attribute) : G(G), hasWeightedEdges(false), hasMSF(false), hasAttribute(false) {
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

#endif // NETWORKIT_GRAPH_RANDOM_MAXIMUM_SPANNING_FOREST_HPP_
