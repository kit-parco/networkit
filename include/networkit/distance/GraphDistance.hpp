// no-networkit-format
/*
 *  GraphDistance.hpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_DISTANCE_GRAPH_DISTANCE_HPP_
#define NETWORKIT_DISTANCE_GRAPH_DISTANCE_HPP_

#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

// TODO: inherit from NodeDistance
/**
 * @ingroup distance
 */
class GraphDistance final {
public:

    /** Default destructor */
    virtual ~GraphDistance() = default;

    /**
     * Returns the distance between @a u and @a v in Graph @a g i.e., the length of the shortest path
     * between the two. Zero if u = v, maximal possible value if no path exists.
     *
     * @param g The graph.
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @return The distance between @a u and @a v.
     */
    edgeweight weightedDistance(const Graph& g, node u, node v) const;

    /**
     * Returns the number of edges on shortest unweighted path between @a u and @a v in Graph @a g.
     *
     * @param g The graph.
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @return The number of edges between @a u and @a v.
     */
    count unweightedDistance(const Graph& g, node u, node v) const;
};

} /* namespace NetworKit */
#endif // NETWORKIT_DISTANCE_GRAPH_DISTANCE_HPP_
