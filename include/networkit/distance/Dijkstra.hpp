// no-networkit-format
/*
 * Dijkstra.hpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning, Christian Staudt
 */

#ifndef NETWORKIT_DISTANCE_DIJKSTRA_HPP_
#define NETWORKIT_DISTANCE_DIJKSTRA_HPP_

#include <tlx/container/d_ary_addressable_int_heap.hpp>

#include <networkit/distance/SSSP.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Dijkstra's SSSP algorithm.
 */
class Dijkstra final : public SSSP {

public:
    /**
     * Creates the Dijkstra class for @a G and the source node @a source.
     *
     * @param G The graph.
     * @param source The source node.
     * @param storePaths Paths are reconstructable and the number of paths is
     *        stored.
     * @param storeNodesSortedByDistance Store a vector of nodes ordered in
     *        increasing distance from the source.
     * @param target The target node.
     */
    Dijkstra(const Graph &G, node source, bool storePaths = true,
             bool storeNodesSortedByDistance = false, node target = none);

    /**
     * Performs the Dijkstra SSSP algorithm on the graph given in the
     * constructor.
     */
    void run() override;

private:
    struct Compare {
    public:
        Compare(const std::vector<double> &dist_) : dist(dist_) {}

        bool operator()(node x, node y) const { return dist[x] < dist[y]; }

    private:
        const std::vector<double> &dist;
    };

    tlx::d_ary_addressable_int_heap<node, 2, Compare> heap;
};

} /* namespace NetworKit */
#endif // NETWORKIT_DISTANCE_DIJKSTRA_HPP_
