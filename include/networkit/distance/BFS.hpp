// no-networkit-format
/*
 * BFS.hpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_DISTANCE_BFS_HPP_
#define NETWORKIT_DISTANCE_BFS_HPP_

#include <networkit/distance/SSSP.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * The BFS class is used to do a breadth-first search on a Graph from a given
 * source node.
 */
class BFS final : public SSSP {

public:
    /**
     * Constructs the BFS class for @a G and source node @a source.
     *
     * @param G The graph
     * @param source The source node of the breadth-first search
     * @param storePaths Paths are reconstructable and the number of paths is
     * stored.
     * @param storeNodesSortedByDistance Store a vector of nodes ordered in
     * increasing distance from the source.
     * @param target The target node.
     */
    BFS(const Graph &G, node source, bool storePaths = true,
        bool storeNodesSortedByDistance = false, node target = none);

    /**
     * Breadth-first search from @a source.
     */
    void run() override;
};
} /* namespace NetworKit */
#endif // NETWORKIT_DISTANCE_BFS_HPP_
