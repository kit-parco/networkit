/*
 * GraphCoarsening.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_COARSENING_GRAPH_COARSENING_HPP_
#define NETWORKIT_COARSENING_GRAPH_COARSENING_HPP_

#include <map>
#include <vector>

#include <networkit/graph/Graph.hpp>
#include <networkit/base/Algorithm.hpp>

namespace NetworKit {

/**
 * @ingroup coarsening
 * Abstract base class for graph coarsening/contraction algorithms.
 */
class GraphCoarsening : public Algorithm {

public:

    GraphCoarsening(const Graph& G);

    virtual ~GraphCoarsening() = default;

    virtual void run() = 0;

    const Graph& getCoarseGraph() const;

    /**
     * Get mapping from fine to coarse node.
     */
    const std::vector<node>& getFineToCoarseNodeMapping() const;

    /**
     * Get mapping from coarse node to collection of fine nodes.
     */
    std::map<node, std::vector<node> > getCoarseToFineNodeMapping() const;

protected:
    const Graph* G;
    Graph Gcoarsened;
    std::vector<node> nodeMapping;

};

} // namespace


#endif // NETWORKIT_COARSENING_GRAPH_COARSENING_HPP_
