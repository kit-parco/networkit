/*
 * DynSSSP.h
 *
 *  Created on: 17.07.2014
 *      Author: cls, ebergamini
 */

#ifndef NETWORKIT_DISTANCE_DYN_SSSP_HPP_
#define NETWORKIT_DISTANCE_DYN_SSSP_HPP_

#include <set>

#include <networkit/graph/Graph.hpp>
#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/distance/SSSP.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Interface for dynamic single-source shortest path algorithms.
 */
class DynSSSP: public SSSP, public DynAlgorithm {

friend class DynApproxBetweenness;
friend class DynApproxBetweennessImpl;
friend class DynApproxBetweennessImplDir;

public:
    /**
     * The algorithm computes a dynamic SSSP starting from the specified
     * source vertex.
     *
     * @param	graph		input graph.
     * @param	source	    source vertex.
     * @param   storePredecessors   keep track of the lists of predecessors?
     */
    DynSSSP(const Graph& G, node source, bool storePredecessors = true, node target = none);

    ~DynSSSP() override = default;

    /**
    * Returns true or false depending on whether the node previoulsy specified
    * with setTargetNode has been modified by the udate or not.
    *
    * @param batch The batch of edge insertions.
    */
    bool modified();
    /**
    * Set a target node to be `observed` during the update. If a node t is set as
    * target before the update, the function modified() will return true or false
    * depending on whether node t has been modified by the update.
    *
    * @param t Node to be `observed`.
    */
    void setTargetNode(node t = 0);

    /**
     * Returns the predecessor nodes of @a t on all shortest paths from source to @a t.
     * @param t Target node.
     * @return The predecessors of @a t on all shortest paths from source to @a t.
     */
    const std::vector<node>& getPredecessors(node t) const;

protected:
    bool storePreds = true;
    bool mod = false;
    node target;
};

inline bool DynSSSP::modified() {
    return mod;
}

inline void DynSSSP::setTargetNode(const node t) {
    target = t;
}

inline const std::vector<node>& DynSSSP::getPredecessors(node t) const {
    if (! storePreds) {
        throw std::runtime_error("predecessors have not been stored");
    }
    return previous[t];
}

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_DYN_SSSP_HPP_
