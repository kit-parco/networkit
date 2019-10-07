/*
 * Graph.h
 *
 *  Created on: 01.06.2014
 *      Author: Christian Staudt (christian.staudt@kit.edu), Klara Reichard
 * (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <algorithm>
#include <functional>
#include <queue>
#include <stack>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/FunctionTraits.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/viz/Point.hpp>
#include <networkit/graph/Coordinates.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

#ifdef NETWORKIT_U32_NODES
	/**
	 * Data type of nodes inside adjacency lists.
	 * This is an internal detail of NetworKit; there is no need for users to use this type.
	 */
	using storednode = uint32_t;
#else
	/**
	 * Data type of nodes inside adjacency lists.
	 * This is an internal detail of NetworKit; there is no need for users to use this type.
	 */
	using storednode = node;
#endif

/**
 * A weighted edge used for the graph constructor with
 * initializer list syntax.
 */
struct WeightedEdge {
    node u, v;
    edgeweight weight;

    WeightedEdge(node u, node v, edgeweight w) : u(u), v(v), weight(w) {}
};
inline bool operator<(const WeightedEdge &e1, const WeightedEdge &e2) {
    return e1.weight < e2.weight;
}
struct Edge {
    node u, v;

    Edge(node _u, node _v, bool sorted = false) {
        if (sorted) {
            u = std::min(_u, _v);
            v = std::max(_u, _v);
        } else {
            u = _u;
            v = _v;
        }
    }
};
inline bool operator==(const Edge &e1, const Edge &e2) {
    return e1.u == e2.u && e1.v == e2.v;
}
struct Unsafe {};
static constexpr Unsafe unsafe {};
} // namespace NetworKit

namespace std {
template <> struct hash<NetworKit::Edge> {
    size_t operator()(const NetworKit::Edge &e) const {
        return hash_node(e.u) ^ hash_node(e.v);
    }

    hash<NetworKit::node> hash_node;
};
} // namespace std

namespace NetworKit {

// forward declaration to randomization/CurveballImpl.h
namespace CurveballDetails {
class CurveballMaterialization;
}

/**
 * @ingroup graph
 * A graph (with optional weights) and parallel iterator methods.
 */
class Graph final {

    friend class ParallelPartitionCoarsening;
    friend class GraphBuilder;
    friend class CurveballDetails::CurveballMaterialization;

  private:
    // graph attributes
    //!< unique graph id, starts at 0
    count id;
    //!< name of the graph, initially G#ID
    std::string name;

    // scalars
    //!< current number of nodes
    count n;
    //!< current number of edges
    count m;

    //!< current number of self loops, edges which have the same origin and
    //!< target
    count storedNumberOfSelfLoops;

    //!< current upper bound of node ids, z will be the id of the next node
    node z;
    //!< current upper bound of edge ids, will be the id of the next edge
    edgeid omega;
    //!< current time step
    count t;

    //!< true if the graph is weighted, false otherwise
    bool weighted;
    //!< true if the graph is directed, false otherwise
    bool directed;
    //!< true if edge ids have been assigned
    bool edgesIndexed;

    // per node data
    //!< exists[v] is true if node v has not been removed from the graph
    std::vector<bool> exists;
    //!< coordinates of nodes (if present)
    Coordinates<float> coordinates;

    //!< only used for directed graphs, inEdges[v] contains all nodes u that
    //!< have an edge (u, v)
    std::vector<std::vector<storednode>> inEdges;
    //!< (outgoing) edges, for each edge (u, v) v is saved in outEdges[u] and
    //!< for undirected also u in outEdges[v]
    std::vector<std::vector<storednode>> outEdges;

    //!< only used for directed graphs, same schema as inEdges
    std::vector<std::vector<edgeweight>> inEdgeWeights;
    //!< same schema (and same order!) as outEdges
    std::vector<std::vector<edgeweight>> outEdgeWeights;

    //!< only used for directed graphs, same schema as inEdges
    std::vector<std::vector<edgeid>> inEdgeIds;
    //!< same schema (and same order!) as outEdges
    std::vector<std::vector<edgeid>> outEdgeIds;

    /**
     * Returns the next unique graph id.
     */
    count getNextGraphId();

    /**
     * Returns the index of node u in the array of incoming edges of node v.
     * (for directed graphs inEdges is searched, while for indirected outEdges
     * is searched, which gives the same result as indexInOutEdgeArray).
     */
    index indexInInEdgeArray(node v, node u) const;

    /**
     * Returns the index of node v in the array of outgoing edges of node u.
     */
    index indexInOutEdgeArray(node u, node v) const;

    /**
     * Computes the weighted in/out degree of a graph.
     *
     * @param inDegree whether to compute the in degree or the out degree.
     */
    edgeweight computeWeightedDegree(const node &v,
                                     const bool inDegree = false) const;

    /**
     * Computes the maximum in/out degree of the graph.
     *
     * @param inDegree wheter to compute the in degree or the out degree.
     */
    count computeMaxDegree(const bool inDegree = false) const;

    /**
     * Computes the maximum in/out weighted degree of the graph
     *
     * @param inDegree whether to compute the in degree or the out degree
     */
    edgeweight computeMaxWeightedDegree(const bool inDegree = false) const;

    /**
     * Returns the edge weight of the outgoing edge of index i in the outgoing
     * edges of node u
     * @param u The node
     * @param i The index
     * @return The weight of the outgoing edge or defaultEdgeWeight if the graph
     * is unweighted
     */
    template <bool hasWeights>
    inline edgeweight getOutEdgeWeight(node u, index i) const;

    /**
     * Returns the edge weight of the incoming edge of index i in the incoming
     * edges of node u
     *
     * @param u The node
     * @param i The index in the incoming edge array
     * @return The weight of the incoming edge
     */
    template <bool hasWeights>
    inline edgeweight getInEdgeWeight(node u, index i) const;

    /**
     * Returns the edge id of the edge of index i in the outgoing edges of node
     * u
     *
     * @param u The node
     * @param i The index in the outgoing edges
     * @return The edge id
     */
    template <bool graphHasEdgeIds>
    inline edgeid getOutEdgeId(node u, index i) const;

    /**
     * Returns the edge id of the edge of index i in the incoming edges of node
     * u
     *
     * @param u The node
     * @param i The index in the incoming edges of u
     * @return The edge id
     */
    template <bool graphHasEdgeIds>
    inline edgeid getInEdgeId(node u, index i) const;

    /**
     * @brief Returns if the edge (u, v) shall be used in the iteration of all
     * edgesIndexed
     *
     * @param u The source node of the edge
     * @param v The target node of the edge
     * @return If the node shall be used, i.e. if v is not none and in the
     * undirected case if u >= v
     */
    template <bool graphIsDirected>
    inline bool useEdgeInIteration(node u, node v) const;

    /**
     * @brief Implementation of the for loop for outgoing edges of u
     *
     * Note: If all (valid) outgoing edges shall be considered, graphIsDirected
     * needs to be set to true
     *
     * @param u The node
     * @param handle The handle that shall be executed for each edge
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds,
              typename L>
    inline void forOutEdgesOfImpl(node u, L handle) const;

    /**
     * @brief Implementation of the for loop for incoming edges of u
     *
     * For undirected graphs, this is the same as forOutEdgesOfImpl but u and v
     * are changed in the handle
     *
     * @param u The node
     * @param handle The handle that shall be executed for each edge
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds,
              typename L>
    inline void forInEdgesOfImpl(node u, L handle) const;

    /**
     * @brief Implementation of the for loop for all edges, @see forEdges
     *
     * @param handle The handle that shall be executed for all edges
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds,
              typename L>
    inline void forEdgeImpl(L handle) const;

    /**
     * @brief Parallel implementation of the for loop for all edges, @see
     * parallelForEdges
     *
     * @param handle The handle that shall be executed for all edges
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds,
              typename L>
    inline void parallelForEdgesImpl(L handle) const;

    /**
     * @brief Summation variant of the parallel for loop for all edges, @see
     * parallelSumForEdges
     *
     * @param handle The handle that shall be executed for all edges
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds,
              typename L>
    inline double parallelSumForEdgesImpl(L handle) const;

    /*
     * In the following definition, Aux::FunctionTraits is used in order to only
     * execute lambda functions with the appropriate parameters. The
     * decltype-return type is used for determining the return type of the
     * lambda (needed for summation) but also determines if the lambda accepts
     * the correct number of parameters. Otherwise the return type declaration
     * fails and the function is excluded from overload resolution. Then there
     * are multiple possible lambdas with three (third parameter id or weight)
     * and two (second parameter can be second node id or edge weight for
     * neighbor iterators). This is checked using Aux::FunctionTraits and
     * std::enable_if. std::enable_if only defines the type member when the
     * given bool is true, this bool comes from std::is_same which compares two
     * types. The function traits give either the parameter type or if it is out
     * of bounds they define type as void.
     */

    /**
     * Triggers a static assert error when no other method is chosen. Because of
     * the use of "..." as arguments, the priority of this method is lower than
     * the priority of the other methods. This method avoids ugly and unreadable
     * template substitution error messages from the other declarations.
     */
    template <class F, void * = (void *)0>
    typename Aux::FunctionTraits<F>::result_type edgeLambda(F &, ...) const {
        // the strange condition is used in order to delay the evaluation of the
        // static assert to the moment when this function is actually used
        static_assert(
            !std::is_same<F, F>::value,
            "Your lambda does not support the required parameters or the "
            "parameters have the wrong type.");
        return std::declval<typename Aux::FunctionTraits<
            F>::result_type>(); // use the correct return type (this won't
                                // compile)
    }

    /**
     * Calls the given function f if its fourth argument is of the type edgeid
     * and third of type edgeweight Note that the decltype check is not enough
     * as edgeweight can be casted to node and we want to assure that .
     */
    template <
        class F,
        typename std::enable_if<
            (Aux::FunctionTraits<F>::arity >= 3) &&
            std::is_same<edgeweight, typename Aux::FunctionTraits<
                                         F>::template arg<2>::type>::value &&
            std::is_same<edgeid, typename Aux::FunctionTraits<F>::template arg<
                                     3>::type>::value>::type * = (void *)0>
    auto edgeLambda(F &f, node u, node v, edgeweight ew, edgeid id) const
        -> decltype(f(u, v, ew, id)) {
        return f(u, v, ew, id);
    }

    /**
     * Calls the given function f if its third argument is of the type edgeid,
     * discards the edge weight Note that the decltype check is not enough as
     * edgeweight can be casted to node.
     */
    template <
        class F,
        typename std::enable_if<
            (Aux::FunctionTraits<F>::arity >= 2) &&
            std::is_same<edgeid, typename Aux::FunctionTraits<F>::template arg<
                                     2>::type>::value &&
            std::is_same<node, typename Aux::FunctionTraits<F>::template arg<
                                   1>::type>::value /* prevent f(v, weight, eid)
                                                     */
            >::type * = (void *)0>
    auto edgeLambda(F &f, node u, node v, edgeweight, edgeid id) const
        -> decltype(f(u, v, id)) {
        return f(u, v, id);
    }

    /**
     * Calls the given function f if its third argument is of type edgeweight,
     * discards the edge id Note that the decltype check is not enough as node
     * can be casted to edgeweight.
     */
    template <class F,
              typename std::enable_if<
                  (Aux::FunctionTraits<F>::arity >= 2) &&
                  std::is_same<edgeweight,
                               typename Aux::FunctionTraits<F>::template arg<
                                   2>::type>::value>::type * = (void *)0>
    auto edgeLambda(F &f, node u, node v, edgeweight ew, edgeid /*id*/) const -> decltype(f(u, v, ew)) {
        return f(u, v, ew);
    }

    /**
     * Calls the given function f if it has only two arguments and the second
     * argument is of type node, discards edge weight and id Note that the
     * decltype check is not enough as edgeweight can be casted to node.
     */
    template <
        class F,
        typename std::enable_if<
            (Aux::FunctionTraits<F>::arity >= 1) &&
            std::is_same<node, typename Aux::FunctionTraits<F>::template arg<
                                   1>::type>::value>::type * = (void *)0>
    auto edgeLambda(F &f, node u, node v, edgeweight /*ew*/, edgeid /*id*/) const -> decltype(f(u, v)) {
        return f(u, v);
    }

    /**
     * Calls the given function f if it has only two arguments and the second
     * argument is of type edgeweight, discards the first node and the edge id
     * Note that the decltype check is not enough as edgeweight can be casted to
     * node.
     */
    template <class F,
              typename std::enable_if<
                  (Aux::FunctionTraits<F>::arity >= 1) &&
                  std::is_same<edgeweight,
                               typename Aux::FunctionTraits<F>::template arg<
                                   1>::type>::value>::type * = (void *)0>
    auto edgeLambda(F &f, node, node v, edgeweight ew, edgeid /*id*/) const -> decltype(f(v, ew)) {
        return f(v, ew);
    }

    /**
     * Calls the given function f if it has only one argument, discards the
     * first node id, the edge weight and the edge id
     */
    template <class F, void * = (void *)0>
    auto edgeLambda(F &f, node, node v, edgeweight, edgeid) const -> decltype(f(v)) {
        return f(v);
    }

    /**
     * Calls the given BFS handle with distance parameter
     */
    template <class F>
    auto callBFSHandle(F &f, node u, count dist) const -> decltype(f(u, dist)) {
        return f(u, dist);
    }

    /**
     * Calls the given BFS handle without distance parameter
     */
    template <class F>
    auto callBFSHandle(F &f, node u, count) const -> decltype(f(u)) {
        return f(u);
    }

public:

    /**
     * Class to iterate over the in/out neighbors of a node.
     */
    class NeighborIterator {

      public:
        // The value type of the neighbors (i.e. nodes). Returned by
        // operator*().
        using value_type = node;

        // Reference to the value_type, required by STL.
        using reference = value_type &;

        // Pointer to the value_type, required by STL.
        using pointer = value_type *;

        // STL iterator category.
        using iterator_category = std::forward_iterator_tag;

        // Signed integer type of the result of subtracting two pointers,
        // required by STL.
        using difference_type = ptrdiff_t;

        // Own type.
        using self = NeighborIterator;

        NeighborIterator(std::vector<storednode>::const_iterator nodesIter) : nIter(nodesIter) {}

        /**
         * @brief WARNING: This contructor is required for Python and should not be used as the iterator is not initialised.
         */
        NeighborIterator() {}

        NeighborIterator operator++() {
            auto prev = *this;
            ++nIter;
            return prev;
        }

        NeighborIterator operator++(int) {
            ++nIter;
            return *this;
        }

        NeighborIterator operator--() {
            auto prev = *this;
            --nIter;
            return prev;
        }

        NeighborIterator operator--(int) {
            --nIter;
            return *this;
        }

        bool operator==(const NeighborIterator &rhs) const {
            return nIter == rhs.nIter;
        }

        bool operator!=(const NeighborIterator &rhs) const {
            return nIter != rhs.nIter;
        }

        node operator*() const { return *nIter; }

      private:
        std::vector<storednode>::const_iterator nIter;
    };

    /**
     * Class to iterate over the in/out neighbors of a node including the edge
     * weights. Values are std::pair<node, edgeweight>.
     */
    class NeighborWeightIterator {

      public:
        // The value type of the neighbors (i.e. nodes). Returned by
        // operator*().
        using value_type = std::pair<node, edgeweight>;

        // Reference to the value_type, required by STL.
        using reference = value_type &;

        // Pointer to the value_type, required by STL.
        using pointer = value_type *;

        // STL iterator category.
        using iterator_category = std::forward_iterator_tag;

        // Signed integer type of the result of subtracting two pointers,
        // required by STL.
        using difference_type = ptrdiff_t;

        // Own type.
        using self = NeighborWeightIterator;

        NeighborWeightIterator(
            std::vector<storednode>::const_iterator nodesIter,
            std::vector<edgeweight>::const_iterator weightIter)
            : nIter(nodesIter), wIter(weightIter) {}

        NeighborWeightIterator operator++() {
            auto prev = *this;
            ++nIter;
            ++wIter;
            return prev;
        }

        NeighborWeightIterator operator++(int) {
            ++nIter;
            ++wIter;
            return *this;
        }

        NeighborWeightIterator operator--() {
            auto prev = *this;
            --nIter;
            --wIter;
            return prev;
        }

        NeighborWeightIterator operator--(int) {
            --nIter;
            --wIter;
            return *this;
        }

        bool operator==(const NeighborWeightIterator &rhs) const {
            return nIter == rhs.nIter && wIter == rhs.wIter;
        }

        bool operator!=(const NeighborWeightIterator &rhs) const {
            return nIter != rhs.nIter || wIter != rhs.wIter;
        }

        const std::pair<node, edgeweight> operator*() const {
            return std::make_pair(*nIter, *wIter);
        }

      private:
        std::vector<storednode>::const_iterator nIter;
        std::vector<edgeweight>::const_iterator wIter;
    };

    /**
     * Wrapper class to iterate over a range of the neighbors of a node within
     * a for loop.
     */
    template <bool InEdges = false> class NeighborRange {
      public:
        NeighborRange(const Graph &G, node u) : G(&G), u(u){
			assert(G.hasNode(u));
        };

        NeighborRange() : G(nullptr) {};

        NeighborIterator begin() const {
            assert(G != nullptr);
            if (InEdges)
                return NeighborIterator(G->inEdges[u].begin());
            return NeighborIterator(G->outEdges[u].begin());
        }

        NeighborIterator end() const {
            assert(G != nullptr);
            if (InEdges)
                return NeighborIterator(G->inEdges[u].end());
            return NeighborIterator(G->outEdges[u].end());
        }

      private:
        const Graph *G;
        node u;
    };

    using OutNeighborRange = NeighborRange<false>;

    using InNeighborRange = NeighborRange<true>;
    /**
     * Wrapper class to iterate over a range of the neighbors of a node
     * including the edge weights within a for loop.
     * Values are std::pair<node, edgeweight>.
     */
    template <bool InEdges = false> class NeighborWeightRange {

      public:
        NeighborWeightRange(const Graph &G, node u) : G(G), u(u){};

        NeighborWeightIterator begin() const {
            if (InEdges)
                return NeighborWeightIterator(G.inEdges[u].begin(),
                                              G.inEdgeWeights[u].begin());
            return NeighborWeightIterator(G.outEdges[u].begin(),
                                          G.outEdgeWeights[u].begin());
        }

        NeighborWeightIterator end() const {
            if (InEdges)
                return NeighborWeightIterator(G.inEdges[u].end(),
                                              G.inEdgeWeights[u].end());
            return NeighborWeightIterator(G.outEdges[u].end(),
                                          G.outEdgeWeights[u].end());
        }

      private:
        const Graph &G;
        const node u;
    };

    /**
     * Create a graph of @a n nodes. The graph has assignable edge weights if @a
     * weighted is set to <code>true</code>. If @a weighted is set to
     * <code>false</code> each edge has edge weight 1.0 and any other weight
     * assignment will be ignored.
     * @param n Number of nodes.
     * @param weighted If set to <code>true</code>, the graph has edge weights.
     * @param directed If set to @c true, the graph will be directed.
     */
    Graph(count n = 0, bool weighted = false, bool directed = false);

    Graph(const Graph &G, bool weighted, bool directed);

    /**
     * Generate a weighted graph from a list of edges. (Useful for small
     * graphs in unit tests that you do not want to read from a file.)
     *
     * @param[in] edges list of weighted edges
     */
    Graph(std::initializer_list<WeightedEdge> edges);

    /**
     * Create a graph as copy of @a other.
     * @param other The graph to copy.
     */
    Graph(const Graph &other) = default;

    /** Default move constructor */
    Graph(Graph &&other) = default;

    /** Default destructor */
    ~Graph() = default;

    /** Default move assignment operator */
    Graph &operator=(Graph &&other) = default;

    /** Default copy assignment operator */
    Graph &operator=(const Graph &other) = default;

	/**
	 * Reserves memory in the node's edge containers for undirected graphs.
	 *
	 * @param u the node memory should be reserved for
	 * @param size the amount of memory to reserve
	 *
	 * This function is thread-safe if called from different
	 * threads on different nodes.
	 */
	void preallocateUndirected(node u, size_t size);

	/**
	 * Reserves memory in the node's edge containers for directed graphs.
	 *
	 * @param u the node memory should be reserved for
	 * @param inSize the amount of memory to reserve for in edges
	 * @param outSize the amount of memory to reserve for out edges
	 *
	 * This function is thread-safe if called from different
	 * threads on different nodes.
	 */
	void preallocateDirected(node u, size_t outSize, size_t inSize);

    /** EDGE IDS **/

    /**
     * Initially assign integer edge identifiers.
     *
     * @param force Force re-indexing of edges even if they have already been
     * indexed
     */
    void indexEdges(bool force = false);

    /**
     * Checks if edges have been indexed
     *
     * @return bool if edges have been indexed
     */
    bool hasEdgeIds() const { return edgesIndexed; }

    /**
     * Get the id of the given edge.
     */
    edgeid edgeId(node u, node v) const;

    /**
     * Get an upper bound for the edge ids in the graph.
     * @return An upper bound for the edge ids.
     */
    index upperEdgeIdBound() const { return omega; }

    /** GRAPH INFORMATION **/

    /**
     * Get the ID of this graph. The ID is a unique unsigned integer given to
     * every graph on construction.
     */
    count getId() const { return id; }

    /**
     * Return the type of the graph.
     * 		Graph: not weighted, undirected
     * 		WeightedGraph: weighted, undirected
     * 		DirectedGraph: not weighted, directed
     * 		WeightedDirectedGraph: weighted, directed
     */
    std::string typ() const;

    /**
     * Try to save some memory by shrinking internal data structures of the
     * graph. Only run this once you finished editing the graph. Otherwise it
     * will cause unnecessary reallocation of memory.
     */
    void shrinkToFit();

    /**
     * Compacts the adjacency arrays by re-using no longer needed slots from
     * deleted edges.
     */
    void compactEdges();

    /**
     * Sorts the adjacency arrays by node id. While the running time is linear
     * this temporarily duplicates the memory.
     */
    void sortEdges();

    /**
     * Set name of graph to @a name.
     * @param name The name.
     */
    void setName(std::string name) { this->name = name; }

    /*
     * Returns the name of the graph.
     * @return The name of the graph.
     */
    std::string getName() const { return name; }

	/**
	 * Set edge count of the graph to edges.
	 * @param edges the edge count of a graph
	 */
	void setEdgeCount(Unsafe, count edges) { m = edges; }

	void setNumberOfSelfLoops(Unsafe, count loops) { storedNumberOfSelfLoops = loops; }
	/**
     * Returns a string representation of the graph.
     * @return A string representation.
     */
    std::string toString() const;

    /* COPYING */

    /*
     * Copies all nodes to a new graph
     * @return graph with the same nodes.
     */
    Graph copyNodes() const;

    /* NODE MODIFIERS */

    /**
     * Add a new node to the graph and return it.
     * @return The new node.
     */
    node addNode();

    /**
     * DEPRECATED: Coordinates should be handled outside the Graph class
     * like general node attributes.
     *
     * Add a new node to the graph with coordinates @a x and @y and return it.
     */
    // TODO: remove method
    // [[deprecated("Deprecated: Node coordinates should be stored externally
    // like any other node attribute")]]
    node addNode(float x, float y);

    /**
     * Remove a node @a v and all incident edges from the graph.
     *
     * Incoming as well as outgoing edges will be removed.
     *
     * @param u Node.
     */
    void removeNode(node v);

    /**
     * Check if node @a v exists in the graph.
     *
     * @param v Node.
     * @return @c true if @a v exists, @c false otherwise.
     */

    bool hasNode(node v) const { return (v < z) && this->exists[v]; }

    /**
     * Restores a previously deleted node @a v with its previous id in the
     * graph.
     *
     * @param v Node.
     *
     */

    void restoreNode(node v);

    // SET OPERATIONS

    /**
     * Appends another graph to this graph as a new subgraph. Performs node
     * id remapping.
     * @param G [description]
     */
    void append(const Graph &G);

    /**
     * Modifies this graph to be the union of it and another graph.
     * Nodes with the same ids are identified with each other.
     * @param G [description]
     */
    void merge(const Graph &G);

    // SUBGRAPHS

    /**
     * Returns an induced subgraph of this graph (including potential edge weights/directions)
     *
     * There a two relevant sets of nodes:
     *  - Nodes are such passed as arguments
     *  - Neighbors are empty by default.
     *      If includeOutNeighbors is set, it includes all out neighbors of Nodes
     *      If includeInNeighbors is set, it includes all in neighbors of Nodes (relevant only for directed graphs)
     *
     * The subgraph contains all nodes in Nodes + Neighbors and all edge which have one end point in Nodes
     * and the other in Nodes or Neighbors.
     */
    Graph subgraphFromNodes(const std::unordered_set<node> &nodes, bool includeOutNeighbors = false, bool includeInNeighbors = false) const;

    /** NODE PROPERTIES **/

    /**
     * Returns the number of outgoing neighbors of @a v.
     *
     * @param v Node.
     * @return The number of outgoing neighbors.
     */
    count degree(node v) const { return outEdges[v].size(); }

    /**
     * Get the number of incoming neighbors of @a v.
     *
     * @param v Node.
     * @return The number of incoming neighbors.
     * @note If the graph is not directed, the outgoing degree is returned.
     */
    count degreeIn(node v) const {
        return directed ? inEdges[v].size() : outEdges[v].size();
    }

    /**
     * Get the number of outgoing neighbors of @a v.
     *
     * @param v Node.
     * @return The number of outgoing neighbors.
     */
    count degreeOut(node v) const { return outEdges[v].size(); }

    /**
     * Returns the maximum out-degree of the graph.
     *
     * @return The maximum out-degree of the graph.
     */
    count maxDegree() const;

    /**
     * Returns the maximum in-degree of the graph.
     *
     * @return The maximum in-degree of the graph.
     */
    count maxDegreeIn() const;

    /**
     * Check whether @a v is isolated, i.e. degree is 0.
     * @param v Node.
     * @return @c true if the node is isolated (= degree is 0)
     */
    bool isIsolated(node v) const {
        if (!exists[v])
            throw std::runtime_error("Error, the node does not exist!");
        return outEdges[v].empty() && (!directed || inEdges[v].empty());
    }

    /**
     * Returns the weighted degree of @a v.
     *
     * @param v Node.
     * @return Weighted degree of @a v.
     * @note For directed graphs this is the sum of weights of all outgoing
     * edges. of @a v.
     */
    edgeweight weightedDegree(const node &v) const;

    /**
     * Returns the maximum weighted degree of the graph.
     *
     * @return Maximum weighted degree of the graph.
     * @note For directed graphs this is the sum of weights of all outgoing
     * edges.
     */
    edgeweight maxWeightedDegree() const;

    /**
     * Returns the maximum weighted in degree of the graph.
     *
     * @return Maximum weighted in degree of the graph.
     * @note For directed graphs this is the sum of weights of all in-going
     * edges.
     */
    edgeweight maxWeightedDegreeIn() const;

    /**
     * Returns the weighted in-degree of @a v.
     *
     * @param v Node.
     * @return Weighted in-degree of @a v.
     * @note For directed graphs this is the sum of weights of all ingoing
     * edges. of @a v.
     */
    edgeweight weightedDegreeIn(const node &v) const;

    /**
     * Returns the volume of the @a v, which is the weighted degree with
     * self-loops counted twice.
     *
     * @param v Node.
     * @return The volume of the @a v.
     */
    edgeweight volume(node v) const;

    /**
     * Returns a random node of the graph.
     * @return A random node.
     */
    node randomNode() const;

    /**
     * Returns a random neighbor of @a u and @c none if degree is zero.
     *
     * @param u Node.
     * @return A random neighbor of @a u.
     */
    node randomNeighbor(node u) const;

    /* EDGE MODIFIERS */

    /**
     * Insert an edge between the nodes @a u and @a v. If the graph is
     * weighted you can optionally set a weight for this edge. The default
     * weight is 1.0. Note: Multi-edges are not supported and will NOT be
     * handled consistently by the graph data structure.
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param weight Optional edge weight.
     */
    void addEdge(node u, node v, edgeweight ew = defaultEdgeWeight);

    /**
    * Insert an edge between the nodes @a u and @a v. Unline the addEdge function, this function does not not add any information to v. If the graph is
    * weighted you can optionally set a weight for this edge. The default
    * weight is 1.0. Note: Multi-edges are not supported and will NOT be
    * handled consistently by the graph data structure.
    * @param u Endpoint of edge.
    * @param v Endpoint of edge.
    * @param weight Optional edge weight.
    * @param ew Optional edge weight.
    * @param index Optional node index.
    */
    void addPartialEdge(Unsafe, node u, node v, edgeweight ew = defaultEdgeWeight, uint64_t index = 0);

    /**
     * Insert an in edge between the nodes @a u and @a v in a directed graph. If the graph is
     * weighted you can optionally set a weight for this edge. The default
     * weight is 1.0. Note: Multi-edges are not supported and will NOT be
     * handled consistently by the graph data structure.
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param ew Optional edge weight.
     * @param index Optional node index.
     */
	void addPartialInEdge(Unsafe, node u, node v, edgeweight ew = defaultEdgeWeight, uint64_t index = 0);

    /**
     * Insert an out edge between the nodes @a u and @a v in a directed graph. If the graph is
     * weighted you can optionally set a weight for this edge. The default
     * weight is 1.0. Note: Multi-edges are not supported and will NOT be
     * handled consistently by the graph data structure.
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param ew Optional edge weight.
     * @param index Optional node index.
     */
	void addPartialOutEdge(Unsafe, node u, node v, edgeweight ew = defaultEdgeWeight, uint64_t index = 0);


    /**
     * Removes the undirected edge {@a u,@a v}.
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     */
    void removeEdge(node u, node v);

    /**
     * Efficiently removes all the edges adjacent to a set of nodes that is
     * not connected to the rest of the graph. This is meant to optimize the
     * Kadabra algorithm.
     * @param nodesInSet vector of nodes that form a connected component that
     * is isolated from the rest of the graph.
     */
    void removeEdgesFromIsolatedSet(const std::vector<node> &nodesInSet);

    /**
     * Removes all the edges in the graph.
     */
    void removeAllEdges();

    /**
     * Removes edges adjacent to a node according to a specific criterion.
     *
     * @param u The node whose adjacent edges shall be removed.
     * @param condition A function that takes a node as an input and returns a
     * bool. If true the edge (u, v) is removed.
     * @param edgesIn Whether in-going or out-going edges shall be removed.
     * @return std::pair<count, count> The number of removed edges (first) and the number of removed
     * self-loops (second).
     */
    template <typename Condition>
    std::pair<count, count> removeAdjacentEdges(node u, Condition condition, bool edgesIn = false);

    /**
     * Removes all self-loops in the graph.
     */
    void removeSelfLoops();

    /**
     * Removes all multi-edges in the graph.
     */
    void removeMultiEdges();

    /**
     * Changes the edges {@a s1, @a t1} into {@a s1, @a t2} and the edge {@a
     * s2,
     * @a t2} into {@a s2, @a t1}.
     *
     * If there are edge weights or edge ids, they are preserved. Note that no
     * check is performed if the swap is actually possible, i.e. does not
     * generate duplicate edges.
     *
     * @param s1 The first source
     * @param t1 The first target
     * @param s2 The second source
     * @param t2 The second target
     */
    void swapEdge(node s1, node t1, node s2, node t2);

    /**
     * Checks if undirected edge {@a u,@a v} exists in the graph.
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @return <code>true</code> if the edge exists, <code>false</code>
     * otherwise.
     */
    bool hasEdge(node u, node v) const;

    /**
     * Returns a random edge. By default a random node u is chosen and then
     * some random neighbor v. So the probability of choosing (u, v) highly
     * depends on the degree of u. Setting uniformDistribution to true, will
     * give you a real uniform distributed edge, but will be slower.
     * Exp. time complexity: O(1) for uniformDistribution = false, O(n) otherwise.
     */
    std::pair<node, node> randomEdge(bool uniformDistribution = false) const;

    /**
     * Returns a vector with nr random edges. The edges are chosen uniform
     * random.
     */
    std::vector<std::pair<node, node>> randomEdges(count nr) const;

    /* GLOBAL PROPERTIES */

    /**
     * Returns <code>true</code> if this graph supports edge weights other
     * than 1.0.
     * @return <code>true</code> if this graph supports edge weights other
     * than 1.0.
     */
    bool isWeighted() const { return weighted; }

    /**
     * Return @c true if this graph supports directed edges.
     * @return @c true if this graph supports directed edges.
     */
    bool isDirected() const { return directed; }

    /**
     * Return <code>true</code> if graph contains no nodes.
     * @return <code>true</code> if graph contains no nodes.
     */
    bool isEmpty() const { return n == 0; }

    /**
     * Return the number of nodes in the graph.
     * @return The number of nodes.
     */
    count numberOfNodes() const { return n; }

    /**
     * Return the number of edges in the graph.
     * @return The number of edges.
     */
    count numberOfEdges() const { return m; }

    /**
     * @return a pair (n, m) where n is the number of nodes and m is the
     * number of edges
     */
    std::pair<count, count> const size() const { return {n, m}; };

    /**
     * @return the density of the graph
     */
    double density() const {
        count n = numberOfNodes();
        count m = numberOfEdges();
        count loops = numberOfSelfLoops();
        m -= loops;
        double d;
        if (isDirected()) {
            d = m / (double)(n * (n - 1));
        } else {
            d = (2 * m) / (double)(n * (n - 1));
        }
        return d;
    }

    /**
     * Return the number of loops {v,v} in the graph.
     * @return The number of loops.
     * @note This involves calculation, so store result if needed multiple
     * times.
     */
    count numberOfSelfLoops() const;

    /**
     * Get an upper bound for the node ids in the graph.
     * @return An upper bound for the node ids.
     */
    index upperNodeIdBound() const { return z; }

    /**
     * Check for invalid graph states, such as multi-edges.
     * @return False if the graph is in invalid state.
     */
    bool checkConsistency() const;

    /* DYNAMICS */

    /**
     * Trigger a time step - increments counter.
     */
    void timeStep() { t++; }

    /**
     * Get time step counter.
     * @return Time step counter.
     */
    count time() { return t; }

    /* COORDINATES */

    /**
     * DEPRECATED: Coordinates should be handled outside the Graph class
     * like general node attributes.
     *
     * Sets the coordinate of @a v to @a value.
     *
     * @param v Node.
     * @param value The coordinate of @a v.
     */
    // TODO: remove method
    // [[deprecated("Deprecated: Node coordinates should be stored externally
    // like any other node attribute")]]
    void setCoordinate(node v, Point<float> value) {
        coordinates.setCoordinate(v, value);
    }

    /**
     * DEPRECATED: Coordinates should be handled outside the Graph class
     * like general node attributes.
     *
     * Get the coordinate of @a v.
     * @param v Node.
     * @return The coordinate of @a v.
     */
    // TODO: remove method
    // [[deprecated("Deprecated: Node coordinates should be stored externally
    // like any other node attribute")]]
    Point<float> &getCoordinate(node v) { return coordinates.getCoordinate(v); }

    /**
     * DEPRECATED: Coordinates should be handled outside the Graph class
     * like general node attributes.
     *
     * Get minimum coordinate of all coordinates with respect to dimension @a
     * dim.
     * @param dim The dimension to search for minimum.
     * @return The minimum coordinate in dimension @a dim.
     */
    // TODO: remove method
    // [[deprecated("Deprecated: Node coordinates should be stored externally
    // like any other node attribute")]]
    float minCoordinate(count dim) { return coordinates.minCoordinate(dim); }

    /**
     * DEPRECATED: Coordinates should be handled outside the Graph class
     * like general node attributes.
     *
     * Get maximum coordinate of all coordinates with respect to dimension @a
     * dim.
     * @param dim The dimension to search for maximum.
     * @return The maximum coordinate in dimension @a dim.
     */
    // TODO: remove method
    // [[deprecated("Deprecated: Node coordinates should be stored externally
    // like any other node attribute")]]
    float maxCoordinate(count dim) { return coordinates.maxCoordinate(dim); }

    /**
     * DEPRECATED: Coordinates should be handled outside the Graph class
     * like general node attributes.
     *
     * Initializes the coordinates for the nodes in graph.
     * @note This has to be called once and before you set coordinates. Call
     * this method again if new nodes have been added.
     */
    // TODO: remove method
    // [[deprecated("Deprecated: Node coordinates should be stored externally
    // like any other node attribute")]]
    void initCoordinates() { coordinates.init(z); }

    /* EDGE ATTRIBUTES */

    /**
     * Return edge weight of edge {@a u,@a v}. Returns 0 if edge does not
     * exist. BEWARE: Running time is \Theta(deg(u))!
     *
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @return Edge weight of edge {@a u,@a v} or 0 if edge does not exist.
     */
    edgeweight weight(node u, node v) const;

    /**
     * Set the weight of an edge. If the edge does not exist,
     * it will be inserted.
     *
     * @param[in]	u	endpoint of edge
     * @param[in]	v	endpoint of edge
     * @param[in]	weight	edge weight
     */
    void setWeight(node u, node v, edgeweight ew);

    /**
     * Increase the weight of an edge. If the edge does not exist,
     * it will be inserted.
     *
     * @param[in]	u	endpoint of edge
     * @param[in]	v	endpoint of edge
     * @param[in]	weight	edge weight
     */
    void increaseWeight(node u, node v, edgeweight ew);

    /* SUMS */

    /**
     * Returns the sum of all edge weights.
     * @return The sum of all edge weights.
     */
    edgeweight totalEdgeWeight() const;

    /* Collections */

    /**
     * Get list of all nodes.
     * @return List of all nodes.
     */
    std::vector<node> nodes() const;

    /**
     * Get list of edges as node pairs.
     * @return List of edges as node pairs.
     */
    std::vector<std::pair<node, node>> edges() const;

    /**
     * Get list of neighbors of @a u.
     *
     * @param u Node.
     * @return List of neighbors of @a u.
     */
     std::vector<node> TLX_DEPRECATED(neighbors(node u) const);

    /**
     * Get an iterable range over the neighbors of @a.
     *
     * @param u Node.
     * @return Iterator range over the neighbors of @a.
     */
    NeighborRange<false> neighborRange(node u) const {
        assert(exists[u]);
        return NeighborRange<false>(*this, u);
    }

    /**
     * Get an iterable range over the neighbors of @a u including the edge
     * weights.
     *
     * @param u Node.
     * @return Iterator range over pairs of neighbors of @a and corresponding
     * edge weights.
     */
    NeighborWeightRange<false> weightNeighborRange(node u) const {
        assert(isWeighted());
        assert(exists[u]);
        return NeighborWeightRange<false>(*this, u);
    }

    /**
     * Get an iterable range over the in-neighbors of @a.
     *
     * @param u Node.
     * @return Iterator range over pairs of in-neighbors of @a.
     */
    NeighborRange<true> inNeighborRange(node u) const {
        assert(isDirected());
        assert(exists[u]);
        return NeighborRange<true>(*this, u);
    }

    /**
     * Get an iterable range over the in-neighbors of @a u including the
     * edge weights.
     *
     * @param u Node.
     * @return Iterator range over pairs of in-neighbors of @a and corresponding
     * edge weights.
     */
    NeighborWeightRange<true> weightInNeighborRange(node u) const {
        assert(isDirected() && isWeighted());
        assert(exists[u]);
        return NeighborWeightRange<true>(*this, u);
    }

    /**
     * Get i-th (outgoing) neighbor of @a u.
     * WARNING: This function is deprecated or only temporary.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return @a i -th (outgoing) neighbor of @a u, or @c none if no such
     * neighbor exists.
     */
    template <bool graphIsDirected> node getIthNeighbor(node u, index i) const {
        node v = outEdges[u][i];
        if (useEdgeInIteration<graphIsDirected>(u, v))
            return v;
        else
            return none;
    }

    /* Derivative Graphs */

    /**
     * Return an undirected version of this graph.
     *
     * @return undirected graph.
     */
    Graph toUndirected() const;

    /**
     * Return an unweighted version of this graph.
     *
     * @return unweighted graph.
     */
    Graph toUnweighted() const;

    /**
     * Return the transpose of this graph. The graph must be directed.
     *
     * @return transpose of the graph.
     */
    Graph transpose() const;

    /* NODE ITERATORS */

    /**
     * Iterate over all nodes of the graph and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L> void forNodes(L handle) const;

    /**
     * Iterate randomly over all nodes of the graph and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L> void parallelForNodes(L handle) const;

    /** Iterate over all nodes of the graph and call @a handle (lambda
     * closure) as long as @a condition remains true. This allows for breaking
     * from a node loop.
     *
     * @param condition Returning <code>false</code> breaks the loop.
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename C, typename L>
    void forNodesWhile(C condition, L handle) const;

    /**
     * Iterate randomly over all nodes of the graph and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L> void forNodesInRandomOrder(L handle) const;

    /**
     * Iterate in parallel over all nodes of the graph and call handler
     * (lambda closure). Using schedule(guided) to remedy load-imbalances due
     * to e.g. unequal degree distribution.
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L> void balancedParallelForNodes(L handle) const;

    /**
     * Iterate over all undirected pairs of nodes and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameters <code>(node, node)</code>.
     */
    template <typename L> void forNodePairs(L handle) const;

    /**
     * Iterate over all undirected pairs of nodes in parallel and call @a
     * handle (lambda closure).
     *
     * @param handle Takes parameters <code>(node, node)</code>.
     */
    template <typename L> void parallelForNodePairs(L handle) const;

    /* EDGE ITERATORS */

    /**
     * Iterate over all edges of the const graph and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameters <code>(node, node)</code>, <code>(node,
     * node, edgweight)</code>, <code>(node, node, edgeid)</code> or
     * <code>(node, node, edgeweight, edgeid)</code>.
     */
    template <typename L> void forEdges(L handle) const;

    /**
     * Iterate in parallel over all edges of the const graph and call @a
     * handle (lambda closure).
     *
     * @param handle Takes parameters <code>(node, node)</code> or
     * <code>(node, node, edgweight)</code>, <code>(node, node, edgeid)</code>
     * or <code>(node, node, edgeweight, edgeid)</code>.
     */
    template <typename L> void parallelForEdges(L handle) const;

    /* NEIGHBORHOOD ITERATORS */

    /**
     * Iterate over all neighbors of a node and call @a handle (lamdba
     * closure).
     *
     * @param u Node.
     * @param handle Takes parameter <code>(node)</code> or <code>(node,
     * edgeweight)</code> which is a neighbor of @a u.
     * @note For directed graphs only outgoing edges from @a u are considered.
     * A node is its own neighbor if there is a self-loop.
     *
     */
    template <typename L> void forNeighborsOf(node u, L handle) const;

    /**
     * Iterate over all incident edges of a node and call @a handle (lamdba
     * closure).
     *
     * @param u Node.
     * @param handle Takes parameters <code>(node, node)</code>, <code>(node,
     * node, edgeweight)</code>, <code>(node, node, edgeid)</code> or
     * <code>(node, node, edgeweight, edgeid)</code> where the first node is
     * @a u and the second is a neighbor of @a u.
     * @note For undirected graphs all edges incident to @a u are also
     * outgoing edges.
     */
    template <typename L> void forEdgesOf(node u, L handle) const;

    /**
     * Iterate over all neighbors of a node and call handler (lamdba closure).
     * For directed graphs only incoming edges from u are considered.
     */
    template <typename L> void forInNeighborsOf(node u, L handle) const;

    /**
     * Iterate over all incoming edges of a node and call handler (lamdba
     * closure).
     * @note For undirected graphs all edges incident to u are also incoming
     * edges.
     *
     * Handle takes parameters (u, v) or (u, v, w) where w is the edge weight.
     */
    template <typename L> void forInEdgesOf(node u, L handle) const;

    /* REDUCTION ITERATORS */

    /**
     * Iterate in parallel over all nodes and sum (reduce +) the values
     * returned by the handler
     */
    template <typename L> double parallelSumForNodes(L handle) const;

    /**
     * Iterate in parallel over all edges and sum (reduce +) the values
     * returned by the handler
     */
    template <typename L> double parallelSumForEdges(L handle) const;

    /* GRAPH SEARCHES */

    /**
     * Iterate over nodes in breadth-first search order starting from r until
     * connected component of r has been visited.
     *
     * @param r Node.
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L> void BFSfrom(node r, L handle) const;
    template <typename L>
    void BFSfrom(const std::vector<node> &startNodes, L handle) const;

    template <typename L> void BFSEdgesFrom(node r, L handle) const;

    /**
     * Iterate over nodes in depth-first search order starting from r until
     * connected component of r has been visited.
     *
     * @param r Node.
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L> void DFSfrom(node r, L handle) const;

    template <typename L> void DFSEdgesFrom(node r, L handle) const;
};

/* NODE ITERATORS */

template <typename L> void Graph::forNodes(L handle) const {
    for (node v = 0; v < z; ++v) {
        if (exists[v]) {
            handle(v);
        }
    }
}

template <typename L> void Graph::parallelForNodes(L handle) const {
#pragma omp parallel for
    for (omp_index v = 0; v < static_cast<omp_index>(z); ++v) {
        if (exists[v]) {
            handle(v);
        }
    }
}

template <typename C, typename L>
void Graph::forNodesWhile(C condition, L handle) const {
    for (node v = 0; v < z; ++v) {
        if (exists[v]) {
            if (!condition()) {
                break;
            }
            handle(v);
        }
    }
}

template <typename L> void Graph::forNodesInRandomOrder(L handle) const {
    std::vector<node> randVec = nodes();
    std::shuffle(randVec.begin(), randVec.end(), Aux::Random::getURNG());
    for (node v : randVec) {
        handle(v);
    }
}

template <typename L> void Graph::balancedParallelForNodes(L handle) const {
// TODO: define min block size (and test it!)
#pragma omp parallel for schedule(guided)
    for (omp_index v = 0; v < static_cast<omp_index>(z); ++v) {
        if (exists[v]) {
            handle(v);
        }
    }
}

template <typename L> void Graph::forNodePairs(L handle) const {
    for (node u = 0; u < z; ++u) {
        if (exists[u]) {
            for (node v = u + 1; v < z; ++v) {
                if (exists[v]) {
                    handle(u, v);
                }
            }
        }
    }
}

template <typename L> void Graph::parallelForNodePairs(L handle) const {
#pragma omp parallel for schedule(guided)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        if (exists[u]) {
            for (node v = u + 1; v < z; ++v) {
                if (exists[v]) {
                    handle(u, v);
                }
            }
        }
    }
}

/* EDGE ITERATORS */

/* HELPERS */

template <typename T>
void erase(node u, index idx, std::vector<std::vector<T>> &vec);
// implementation for weighted == true
template <bool hasWeights>
inline edgeweight Graph::getOutEdgeWeight(node u, index i) const {
    return outEdgeWeights[u][i];
}

// implementation for weighted == false
template <>
inline edgeweight Graph::getOutEdgeWeight<false>(node, index) const {
    return defaultEdgeWeight;
}

// implementation for weighted == true
template <bool hasWeights>
inline edgeweight Graph::getInEdgeWeight(node u, index i) const {
    return inEdgeWeights[u][i];
}

// implementation for weighted == false
template <> inline edgeweight Graph::getInEdgeWeight<false>(node, index) const {
    return defaultEdgeWeight;
}

// implementation for hasEdgeIds == true
template <bool graphHasEdgeIds>
inline edgeid Graph::getOutEdgeId(node u, index i) const {
    return outEdgeIds[u][i];
}

// implementation for hasEdgeIds == false
template <> inline edgeid Graph::getOutEdgeId<false>(node, index) const {
    return 0;
}

// implementation for hasEdgeIds == true
template <bool graphHasEdgeIds>
inline edgeid Graph::getInEdgeId(node u, index i) const {
    return inEdgeIds[u][i];
}

// implementation for hasEdgeIds == false
template <> inline edgeid Graph::getInEdgeId<false>(node, index) const {
    return 0;
}

// implementation for graphIsDirected == true
template <bool graphIsDirected>
inline bool Graph::useEdgeInIteration(node /* u */, node v) const {
    return v != none;
}

// implementation for graphIsDirected == false
template <> inline bool Graph::useEdgeInIteration<false>(node u, node v) const {
    return u >= v;
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds,
          typename L>
inline void Graph::forOutEdgesOfImpl(node u, L handle) const {
    for (index i = 0; i < outEdges[u].size(); ++i) {
        node v = outEdges[u][i];

        if (useEdgeInIteration<graphIsDirected>(u, v)) {
            edgeLambda<L>(handle, u, v, getOutEdgeWeight<hasWeights>(u, i),
                          getOutEdgeId<graphHasEdgeIds>(u, i));
        }
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds,
          typename L>
inline void Graph::forInEdgesOfImpl(node u, L handle) const {
    if (graphIsDirected) {
        for (index i = 0; i < inEdges[u].size(); i++) {
            node v = inEdges[u][i];

            if (useEdgeInIteration<true>(u, v)) {
                edgeLambda<L>(handle, u, v, getInEdgeWeight<hasWeights>(u, i),
                              getInEdgeId<graphHasEdgeIds>(u, i));
            }
        }
    } else {
        for (index i = 0; i < outEdges[u].size(); ++i) {
            node v = outEdges[u][i];

            if (useEdgeInIteration<true>(u, v)) {
                edgeLambda<L>(handle, u, v, getOutEdgeWeight<hasWeights>(u, i),
                              getOutEdgeId<graphHasEdgeIds>(u, i));
            }
        }
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds,
          typename L>
inline void Graph::forEdgeImpl(L handle) const {
    for (node u = 0; u < z; ++u) {
        forOutEdgesOfImpl<graphIsDirected, hasWeights, graphHasEdgeIds, L>(
            u, handle);
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds,
          typename L>
inline void Graph::parallelForEdgesImpl(L handle) const {
#pragma omp parallel for schedule(guided)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        forOutEdgesOfImpl<graphIsDirected, hasWeights, graphHasEdgeIds, L>(
            u, handle);
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds,
          typename L>
inline double Graph::parallelSumForEdgesImpl(L handle) const {
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        for (index i = 0; i < outEdges[u].size(); ++i) {
            node v = outEdges[u][i];

            // undirected, do not iterate over edges twice
            // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
            if (useEdgeInIteration<graphIsDirected>(u, v)) {
                sum += edgeLambda<L>(handle, u, v,
                                     getOutEdgeWeight<hasWeights>(u, i),
                                     getOutEdgeId<graphHasEdgeIds>(u, i));
            }
        }
    }

    return sum;
}

template <typename L> void Graph::forEdges(L handle) const {
    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edgeIds
        forEdgeImpl<false, false, false, L>(handle);
        break;

    case 1: // weighted,   undirected, no edgeIds
        forEdgeImpl<false, true, false, L>(handle);
        break;

    case 2: // unweighted, directed, no edgeIds
        forEdgeImpl<true, false, false, L>(handle);
        break;

    case 3: // weighted, directed, no edgeIds
        forEdgeImpl<true, true, false, L>(handle);
        break;

    case 4: // unweighted, undirected, with edgeIds
        forEdgeImpl<false, false, true, L>(handle);
        break;

    case 5: // weighted,   undirected, with edgeIds
        forEdgeImpl<false, true, true, L>(handle);
        break;

    case 6: // unweighted, directed, with edgeIds
        forEdgeImpl<true, false, true, L>(handle);
        break;

    case 7: // weighted,   directed, with edgeIds
        forEdgeImpl<true, true, true, L>(handle);
        break;
    }
}

template <typename L> void Graph::parallelForEdges(L handle) const {
    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edgeIds
        parallelForEdgesImpl<false, false, false, L>(handle);
        break;

    case 1: // weighted,   undirected, no edgeIds
        parallelForEdgesImpl<false, true, false, L>(handle);
        break;

    case 2: // unweighted, directed, no edgeIds
        parallelForEdgesImpl<true, false, false, L>(handle);
        break;

    case 3: // weighted, directed, no edgeIds
        parallelForEdgesImpl<true, true, false, L>(handle);
        break;

    case 4: // unweighted, undirected, with edgeIds
        parallelForEdgesImpl<false, false, true, L>(handle);
        break;

    case 5: // weighted,   undirected, with edgeIds
        parallelForEdgesImpl<false, true, true, L>(handle);
        break;

    case 6: // unweighted, directed, with edgeIds
        parallelForEdgesImpl<true, false, true, L>(handle);
        break;

    case 7: // weighted,   directed, with edgeIds
        parallelForEdgesImpl<true, true, true, L>(handle);
        break;
    }
}

/* NEIGHBORHOOD ITERATORS */

template <typename L> void Graph::forNeighborsOf(node u, L handle) const {
    forEdgesOf(u, handle);
}

template <typename L> void Graph::forEdgesOf(node u, L handle) const {
    switch (weighted + 2 * edgesIndexed) {
    case 0: // not weighted, no edge ids
        forOutEdgesOfImpl<true, false, false, L>(u, handle);
        break;

    case 1: // weighted, no edge ids
        forOutEdgesOfImpl<true, true, false, L>(u, handle);
        break;

    case 2: // not weighted, with edge ids
        forOutEdgesOfImpl<true, false, true, L>(u, handle);
        break;

    case 3: // weighted, with edge ids
        forOutEdgesOfImpl<true, true, true, L>(u, handle);
        break;
    }
}

template <typename L> void Graph::forInNeighborsOf(node u, L handle) const {
    forInEdgesOf(u, handle);
}

template <typename L> void Graph::forInEdgesOf(node u, L handle) const {
    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edge ids
        forInEdgesOfImpl<false, false, false, L>(u, handle);
        break;

    case 1: // weighted, undirected, no edge ids
        forInEdgesOfImpl<false, true, false, L>(u, handle);
        break;

    case 2: // unweighted, directed, no edge ids
        forInEdgesOfImpl<true, false, false, L>(u, handle);
        break;

    case 3: // weighted, directed, no edge ids
        forInEdgesOfImpl<true, true, false, L>(u, handle);
        break;

    case 4: // unweighted, undirected, with edge ids
        forInEdgesOfImpl<false, false, true, L>(u, handle);
        break;

    case 5: // weighted, undirected, with edge ids
        forInEdgesOfImpl<false, true, true, L>(u, handle);
        break;

    case 6: // unweighted, directed, with edge ids
        forInEdgesOfImpl<true, false, true, L>(u, handle);
        break;

    case 7: // weighted, directed, with edge ids
        forInEdgesOfImpl<true, true, true, L>(u, handle);
        break;
    }
}

/* REDUCTION ITERATORS */

template <typename L> double Graph::parallelSumForNodes(L handle) const {
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (omp_index v = 0; v < static_cast<omp_index>(z); ++v) {
        if (exists[v]) {
            sum += handle(v);
        }
    }

    return sum;
}

template <typename L> double Graph::parallelSumForEdges(L handle) const {
    double sum = 0.0;

    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edge ids
        sum = parallelSumForEdgesImpl<false, false, false, L>(handle);
        break;

    case 1: // weighted,   undirected, no edge ids
        sum = parallelSumForEdgesImpl<false, true, false, L>(handle);
        break;

    case 2: // unweighted, directed, no edge ids
        sum = parallelSumForEdgesImpl<true, false, false, L>(handle);
        break;

    case 3: // weighted,   directed, no edge ids
        sum = parallelSumForEdgesImpl<true, true, false, L>(handle);
        break;

    case 4: // unweighted, undirected, with edge ids
        sum = parallelSumForEdgesImpl<false, false, true, L>(handle);
        break;

    case 5: // weighted,   undirected, with edge ids
        sum = parallelSumForEdgesImpl<false, true, true, L>(handle);
        break;

    case 6: // unweighted, directed, with edge ids
        sum = parallelSumForEdgesImpl<true, false, true, L>(handle);
        break;

    case 7: // weighted,   directed, with edge ids
        sum = parallelSumForEdgesImpl<true, true, true, L>(handle);
        break;
    }

    return sum;
}

/* GRAPH SEARCHES */

template <typename L> void Graph::BFSfrom(node r, L handle) const {
    std::vector<node> startNodes(1, r);
    BFSfrom(startNodes, handle);
}

template <typename L>
void Graph::BFSfrom(const std::vector<node> &startNodes, L handle) const {
    std::vector<bool> marked(z);
    std::queue<node> q, qNext;
    count dist = 0;
    // enqueue start nodes
    for (node u : startNodes) {
        q.push(u);
        marked[u] = true;
    }
    do {
        node u = q.front();
        q.pop();
        // apply function
        callBFSHandle(handle, u, dist);
        forNeighborsOf(u, [&](node v) {
            if (!marked[v]) {
                qNext.push(v);
                marked[v] = true;
            }
        });
        if (q.empty() && !qNext.empty()) {
            q.swap(qNext);
            ++dist;
        }
    } while (!q.empty());
}

template <typename L> void Graph::BFSEdgesFrom(node r, L handle) const {
    std::vector<bool> marked(z);
    std::queue<node> q;
    q.push(r); // enqueue root
    marked[r] = true;
    do {
        node u = q.front();
        q.pop();
        // apply function
        forNeighborsOf(u, [&](node, node v, edgeweight w, edgeid eid) {
            if (!marked[v]) {
                handle(u, v, w, eid);
                q.push(v);
                marked[v] = true;
            }
        });
    } while (!q.empty());
}

template <typename L> void Graph::DFSfrom(node r, L handle) const {
    std::vector<bool> marked(z);
    std::stack<node> s;
    s.push(r); // enqueue root
    marked[r] = true;
    do {
        node u = s.top();
        s.pop();
        // apply function
        handle(u);
        forNeighborsOf(u, [&](node v) {
            if (!marked[v]) {
                s.push(v);
                marked[v] = true;
            }
        });
    } while (!s.empty());
}

template <typename L> void Graph::DFSEdgesFrom(node r, L handle) const {
    std::vector<bool> marked(z);
    std::stack<node> s;
    s.push(r); // enqueue root
    marked[r] = true;
    do {
        node u = s.top();
        s.pop();
        // apply function
        forNeighborsOf(u, [&](node v) {
            if (!marked[v]) {
                handle(u, v);
                s.push(v);
                marked[v] = true;
            }
        });
    } while (!s.empty());
}

/* EDGE MODIFIERS */

template <typename Condition>
std::pair<count, count> Graph::removeAdjacentEdges(node u, Condition condition, bool edgesIn) {
	count removedEdges = 0;
	count removedSelfLoops = 0;

	// For directed graphs, this function is supposed to be called twice: one to remove out-edges,
	// and one to remove in-edges.
	auto &edges_ = edgesIn ? inEdges[u] : outEdges[u];
	for (index vi = 0; vi < edges_.size();) {
		if (condition(edges_[vi])) {
            const auto isSelfLoop = (edges_[vi] == u);
            removedSelfLoops += isSelfLoop;
            removedEdges += !isSelfLoop;
			edges_[vi] = edges_.back();
			edges_.pop_back();
			if (isWeighted()) {
				auto &weights_ = edgesIn ? inEdgeWeights[u] : outEdgeWeights[u];
				weights_[vi] = weights_.back();
				weights_.pop_back();
			}
			if (hasEdgeIds()) {
				auto &edgeIds_ = edgesIn ? inEdgeIds[u] : outEdgeIds[u];
				edgeIds_[vi] = edgeIds_.back();
				edgeIds_.pop_back();
			}
		} else {
			++vi;
		}
	}

	return {removedEdges, removedSelfLoops};
}


} /* namespace NetworKit */

#endif /* GRAPH_H_ */
