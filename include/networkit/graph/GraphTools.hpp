
#ifndef NETWORKIT_GRAPH_GRAPH_TOOLS_HPP_
#define NETWORKIT_GRAPH_GRAPH_TOOLS_HPP_

#include <unordered_map>
#include <unordered_set>

#include <tlx/define/deprecated.hpp>
#include <tlx/unused.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace GraphTools {

/**
 * Returns the maximum out-degree of the graph.
 *
 * @param G The input graph.
 * @return The maximum out-degree of the graph.
 */
count maxDegree(const Graph &G);

/**
 * Returns the maximum in-degree of the graph.
 *
 * @param G The input graph.
 * @return The maximum in-degree of the graph.
 */
count maxInDegree(const Graph &G);

/**
 * Returns the maximum weighted out-degree of the graph.
 *
 * @param G The input graph.
 * @return Maximum weighted degree of the graph.
 */
edgeweight maxWeightedDegree(const Graph &G);

/**
 * Returns the maximum weighted in-degree of the graph.
 *
 * @param G The input graph.
 * @return Maximum weighted in degree of the graph.
 */
edgeweight maxWeightedInDegree(const Graph &G);

/**
 * Returns a random node of the input graph.
 *
 * @param G The input graph.
 * @return A random node.
 */
node randomNode(const Graph &G);

/**
 * Returns n distinct random nodes of the input graph.
 *
 * @param G The input graph.
 * @param n The number of desired nodes.
 * @return A vector of n disntinct random nodes.
 */
std::vector<node> randomNodes(const Graph &G, count n);

/**
 * Returns a random neighbor of node @a u. Returns none if degree is zero.
 *
 * @param G The input graph.
 * @param u Node.
 * @return A random neighbor of @a u.
 */
node randomNeighbor(const Graph &G, node u);

/**
 * Returns a random edge. By default a random node u is chosen and then
 * some random neighbor v. So the probability of choosing (u, v) highly
 * depends on the degree of u. Setting uniformDistribution to true, will
 * give you a real uniform distributed edge, but will be slower.
 * Exp. time complexity: O(1) for uniformDistribution = false, O(n) otherwise.
 *
 * @param Graph G The input graph.
 * @param bool uniformDistribution Whether the random edge should be extracted uniformly at
 * random.
 * @return std::pair<node, node> A random edge.
 */
std::pair<node, node> randomEdge(const Graph &G, bool uniformDistribution = false);

/**
 * Returns a vector with @a nr random edges. The edges are chosen uniformly
 * random.
 *
 * @param G The input graph.
 * @param nr The number of random edges to be returned.
 * @return std::vector<std::pair<node, node>> Vector with random edges.
 */
std::vector<std::pair<node, node>> randomEdges(const Graph &G, count nr);

/**
 * Efficiently removes all the edges adjacent to a set of nodes that is
 * not connected to the rest of the graph. This is meant to optimize the
 * Kadabra algorithm.
 *
 * @param G The input graph.
 * @param first Start of the range that contains the nodes in the set.
 * @param last End of the range that contains the nodes in the set.
 * is isolated from the rest of the graph.
 */
template <class InputIt>
void removeEdgesFromIsolatedSet(Graph &G, InputIt first, InputIt last) {
    count removedEdges = 0;
    while (first != last) {
        const auto u = *first++;
        removedEdges += G.degree(u);
        G.removePartialOutEdges(unsafe, u);
        if (G.isDirected()) {
            G.removePartialInEdges(unsafe, u);
        }
    }

    G.setEdgeCount(unsafe, G.numberOfEdges() - (G.isDirected() ? removedEdges : removedEdges / 2));
}

/**
 * Returns the number of nodes and the number of edges of the input graph.
 *
 * @param G The input graph.
 * @return std::pair<count, count> with the number of nodes and the number
 * of edges of the input graph.
 */
std::pair<node, node> size(const Graph &G) noexcept;

/**
 * Return the density of the input graph.
 *
 * @param G The input graph.
 *
 * @return double The density of the input graph.
 */
double density(const Graph &G) noexcept;

/**
 * Returns the volume of the input graph.
 *
 * @param G The input graph.
 * @return Volume of the graph.
 */
double volume(const Graph &G);

/**
 * Returns the volume (sum of the out-degree of all nodes) of a range
 * of nodes.
 *
 * @param G The input graph.
 * @param first,last The range of nodes to insert
 * @return Volume of the graph.
 */
template <class InputIt>
double volume(const Graph &G, InputIt first, InputIt last) {
    double outVolume = 0.0;

    // Compute volume of subgraph
    for (; first != last; ++first)
        outVolume += G.weightedDegree(*first, !G.isDirected());

    return outVolume;
}

/**
 * Returns the inVolume (sum of the in-degree of all nodes) of a range
 * of nodes.
 *
 * @param G The input graph.
 * @param first,last The range of nodes to insert
 * @return inVolume of the subgraph.
 */
template <class InputIt>
double inVolume(const Graph &G, InputIt first, InputIt last) {
    double inVolume = 0.0;

    // Compute volume of subgraph
    for (; first != last; ++first)
        inVolume += G.weightedDegreeIn(*first, !G.isDirected());

    return inVolume;
}

/**
 * Copies all nodes of the input graph to a new graph (edges are not copied).
 *
 * @param G The input graph.
 *
 * @return Graph with the same nodes as the input graph (and without any edge).
 */
Graph copyNodes(const Graph &G);

/**
 * Returns an induced subgraph of the input graph (including potential edge weights/directions)
 *
 * The subgraph contains all given nodes and all edges which have both end points in nodes.
 *
 * @param G The input graph.
 * @param nodes The nodes of the induced subgraph.
 *
 * @return Induced subgraph.
 */
Graph subgraphFromNodes(const Graph &G, const std::unordered_set<node> &nodes);

/**
 * Returns an induced subgraph of the input graph (including potential edge weights/directions)
 *
 * @deprecated Use subgraphFromNodesAndEdgesToTheirNeighbors instead.
 *
 * There a two relevant sets of nodes:
 *  - Nodes are such passed as arguments.
 *  - Neighbors are empty by default.
 *
 * The subgraph contains all nodes in Nodes + Neighbors and all edges which have one end point
 * in Nodes and the other in Nodes or Neighbors.
 *
 *
 * @param G The input graph.
 * @param nodes The nodes of the induced subgraph.
 * @param includeOutNeighbors If set to true, all out-neighbors will also be included.
 * @param includeInNeighbors If set to true, all in-neighbors will also be included
 * (relevant only for directed graphs).
 *
 *
 * @return Induced subgraph.
 */
Graph TLX_DEPRECATED(subgraphFromNodes(const Graph &G, const std::unordered_set<node> &nodes,
                                       bool includeOutNeighbors, bool includeInNeighbors = false));

/**
 * Returns an induced subgraph of the input graph (including potential edge weights/directions)
 *
 * The subgraph contains all nodes in the given node range and all edges which have both end points
 * in nodes.
 *
 * @param G The input graph.
 * @param first,last The range of nodes of the induced subgraph.
 * @param compact If the resulting graph shall have compact, continuous node ids, alternatively,
 * node ids of the input graph are kept.
 *
 * @return Induced subgraph.
 */
template <class InputIt>
Graph subgraphFromNodes(const Graph &G, InputIt first, InputIt last, bool compact = false) {
    count subgraphIdBound = 0;
    std::unordered_map<node, node> reverseMapping;

    if (compact) {
        for (InputIt it = first; it != last; ++it) {
            reverseMapping[*it] = subgraphIdBound;
            ++subgraphIdBound;
        }
    } else {
        subgraphIdBound = G.upperNodeIdBound();
    }

    Graph S(subgraphIdBound, G.isWeighted(), G.isDirected());

    if (compact) {
        for (auto nodeIt : reverseMapping) {
            node u = nodeIt.first;
            node localU = nodeIt.second;
            G.forNeighborsOf(u, [&](node v, edgeweight weight) {
                if (!G.isDirected() && u > v)
                    return;

                auto vMapping = reverseMapping.find(v);
                if (vMapping != reverseMapping.end())
                    S.addEdge(localU, vMapping->second, weight);
            });
        }
    } else {
        // First, delete all nodes
        for (node u = 0; u < G.upperNodeIdBound(); ++u) {
            S.removeNode(u);
        }

        // Restore all given nodes
        for (InputIt it = first; it != last; ++it) {
            S.restoreNode(*it);
        }

        G.forEdges([&](node u, node v, edgeweight w) {
            // only include edges if at least one endpoint is in nodes
            if (S.hasNode(u) && S.hasNode(v)) {
                S.addEdge(u, v, w);
            }
        });
    }

    return S;
}

/**
 * Returns an induced subgraph of this graph (including potential edge weights/directions)
 *
 * There a two relevant sets of nodes:
 *  - Nodes are such passed as arguments.
 *  - Neighbors are empty by default.
 *
 * The subgraph contains all nodes in Nodes + Neighbors and all edges which have one end point
 * in Nodes and the other in Nodes or Neighbors.
 *
 * @param G The input graph.
 * @param nodes Nodes of the induced subgraph.
 * @param includeOutNeighbors If set to true, all out-neighbors will also be included.
 * @param includeInNeighbors If set to true, all in-neighbors will also be included
 * (relevant only for directed graphs).
 *
 * @return Induced subgraph.
 */
Graph subgraphAndNeighborsFromNodes(const Graph &G, const std::unordered_set<node> &nodes,
                                    bool includeOutNeighbors = false,
                                    bool includeInNeighbors = false);

/**
 * Returns an undirected copy of the input graph.
 *
 * @param G The input graph.
 *
 * @return Undirected copy of the input graph.
 */
Graph toUndirected(const Graph &G);

/**
 * Return an unweighted copy of the input graph.
 *
 * @param G The input graph.
 *
 * @return Unweighted copy of the input graph.
 */
Graph toUnweighted(const Graph &G);

/**
 * Return a weighted copy of the input graph.
 *
 * @param G The input graph.
 *
 * @return Weighted copy of the input graph.
 */
Graph toWeighted(const Graph &G);

/**
 * Returns the transpose of the input graph. The graph must be directed.
 *
 * @param G The input graph.
 *
 * @return Transpose of the input graph.
 */
Graph transpose(const Graph &G);

/**
 * Appends graph @a G1 to graph @a G as a new subgraph. Performs node id remapping.
 *
 * @param G Graph where @G1 will be appended to.
 * @param G1 Graph that will be appended to @a G.
 */
void append(Graph &G, const Graph &G1);

/**
 * Modifies graph @a G to be the union of it and graph @a G1.
 * Nodes with the same ids are identified with each other.
 *
 * @param G Result of the merge.
 * @param G1 Graph that will be merged with @a G.
 */
void merge(Graph &G, const Graph &G1);

/**
 * Computes a graph with the same structure but with continuous node ids.
 * @param  graph     The graph to be compacted.
 * @param  nodeIdMap The map providing the information about the node ids.
 * @return           Returns a compacted Graph.
 */
Graph getCompactedGraph(const Graph &graph, const std::unordered_map<node, node> &nodeIdMap);

/**
 * Computes a map of node ids.
 * @param	graph	The graph of which the node id map is wanted.
 * @return			Returns the node id map.
 */
std::unordered_map<node, node> getContinuousNodeIds(const Graph &graph);

/**
 * Computes a map of random node ids
 * @param	graph	The graph of which the node id map is wanted.
 * @return		Returns the node id map.
 */
std::unordered_map<node, node> getRandomContinuousNodeIds(const Graph &graph);

/**
 * Inverts a given mapping of node ids from a graph with deleted nodes to continuous node ids.
 * @param 	nodeIdMap	The mapping from node ids with gaps to continuous node ids (i.e. from
 * @getContinuousNodeIds)
 * @param 	G 			The compacted graph (currently only needed for the upper node id bound)
 * @return 				A vector of nodes id where the index is the node id of the compacted graph
 * and the value is the node id of the noncontinuous graph.
 */
std::vector<node> invertContinuousNodeIds(const std::unordered_map<node, node> &nodeIdMap,
                                          const Graph &G);

/**
 * Constructs a new graph that has the same node ids as before it was compacted.
 * @param  invertedIdMap The node id mapping from continuous node ids to noncontinuous node ids.
 * @param  G             The compacted graph.
 * @return               The original graph.
 */
Graph restoreGraph(const std::vector<node> &invertedIdMap, const Graph &G);

/**
 * Sorts the adjacency arrays by increasing or decreasing edge weight. Edge ids are used
 * to break ties.
 *
 * @param G The input graph.
 * @param decreasing If true the adjacency arrays are sorted by non-increasing edge weights, if
 * false, the adjacency arrays are sorted by non-decreasing edge weights. Ties are broken by
 * using node ids.
 */
void sortEdgesByWeight(Graph &G, bool decreasing = false);

/**
 * Rename nodes in a graph using a callback which translates each old id to a new one.
 * For each node u in input graph, oldIdToNew(u) < numNodes.
 *
 * @param graph Input graph.
 * @param numNodes    Number of nodes in the output graph.
 * @param oldIdToNew  Translate old id to new ones. Must be thread-safe
 * @param skipNode    Skip all nodes (and incident edges) for old node
 *                    ids u where deleteNode(u) == true, Must be thread-safe
 * @param preallocate Preallocates memory before adding neighbors
 *                    (Preallocation does not account for deleted nodes
 *                    and hence may need more memory)
 *
 * @node preallocate is currently not implemented
 */
template <typename UnaryIdMapper, typename SkipEdgePredicate>
Graph getRemappedGraph(const Graph &graph, count numNodes, UnaryIdMapper &&oldIdToNew,
                       SkipEdgePredicate &&skipNode, bool preallocate = true) {
    tlx::unused(preallocate); // TODO: Add perallocate as soon as Graph supports it

#ifndef NDEBUG
    graph.forNodes([&](node u) { assert(skipNode(u) || oldIdToNew(u) < numNodes); });
#endif // NDEBUG

    const auto directed = graph.isDirected();
    Graph Gnew(numNodes, graph.isWeighted(), directed);

    graph.forNodes([&](const node u) { // TODO: Make parallel when graph support addHalfEdge
        if (skipNode(u))
            return;

        const node mapped_u = oldIdToNew(u);
        graph.forNeighborsOf(u, [&](node, node v, edgeweight ew) {
            if (!directed && v < u)
                return;
            if (skipNode(v))
                return;

            const node mapped_v = oldIdToNew(v);
            Gnew.addEdge(mapped_u, mapped_v, ew);
        });
    });

    return Gnew;
}

template <typename UnaryIdMapper>
Graph getRemappedGraph(const Graph &graph, count numNodes, UnaryIdMapper &&oldIdToNew,
                       bool preallocate = true) {
    return getRemappedGraph(
        graph, numNodes, std::forward<UnaryIdMapper>(oldIdToNew), [](node) { return false; },
        preallocate);
}

} // namespace GraphTools

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_GRAPH_TOOLS_HPP_
