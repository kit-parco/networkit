// no-networkit-format
/*
 * DynBFS.hpp
 *
 *  Created on: 17.07.2014
 *      Author: cls, ebergamini
 */

#ifndef NETWORKIT_DISTANCE_DYN_BFS_HPP_
#define NETWORKIT_DISTANCE_DYN_BFS_HPP_

#include <networkit/distance/DynSSSP.hpp>


namespace NetworKit {

/**
 * @ingroup distance
 * Dynamic breadth-first search.
 */
class DynBFS final : public DynSSSP {

public:

    /**
     * Creates the object for @a G and source @a s.
     *
     * @param G The graph.
     * @param s The source node.
     * @param storePredecessors keep track of the lists of predecessors?
     */
    DynBFS(const Graph& G, node s, bool storePredecessors = true);

    void run() override;

    /** Updates the distances after an edge insertion.*/
    void update(GraphEvent e) override;

    /** Updates the distances after a batch of edge insertions.*/
    void updateBatch(const std::vector<GraphEvent>& batch) override;

    /* Returns the number of shortest paths to node t.*/
    bigfloat getNumberOfPaths(node t) const;

private:
    enum Color {WHITE, BLACK, GRAY};
    std::vector<Color> color;
    count maxDistance;

};

inline bigfloat DynBFS::getNumberOfPaths(node t) const {
    return npaths[t];
}

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_DYN_BFS_HPP_
