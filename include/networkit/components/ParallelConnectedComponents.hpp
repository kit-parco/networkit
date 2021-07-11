// no-networkit-format
/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#ifndef NETWORKIT_COMPONENTS_PARALLEL_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_PARALLEL_CONNECTED_COMPONENTS_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/base/Algorithm.hpp>

namespace NetworKit {

/**
 * @ingroup components
 * Determines the connected components of an undirected graph.
 */
class ParallelConnectedComponents final : public Algorithm {
public:

    /**
     * @param[in] G Graph for which connected components shall be computed.
     * @param[in] coarsening Specifies whether the main algorithm based on
     *   label propagation (LP) shall work recursively (true) or not (false) by
     *   coarsening/contracting an LP-computed clustering. Defaults to true
     *   since we saw positive effects in terms of running time for many networks.
     *   Beware of possible memory implications.
     */
    ParallelConnectedComponents(const Graph& G, bool coarsening = true);

    /**
     * This method determines the connected components for the graph g.
     */
    void runSequential();

    /**
     * This method determines the connected components for the graph g.
     */
    void run() override;

    /**
     * This method returns the number of connected components.
     */
    count numberOfComponents() const;

    /**
     * This method returns the the component in which node query is situated.
     *
     * @param[in]	query	the node whose component is asked for
     */
    count componentOfNode(node u) const;


    /**
     * Return a Partition that represents the components
     */
    Partition getPartition() const;

    /**
     * @return Vector of components, each stored as (unordered) set of nodes.
     */
    std::vector<std::vector<node>> getComponents() const;

private:
    const Graph* G;
    Partition component;
    bool coarsening;
};

}


#endif // NETWORKIT_COMPONENTS_PARALLEL_CONNECTED_COMPONENTS_HPP_
