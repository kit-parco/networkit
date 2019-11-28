/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#ifndef NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_HPP_

#include <cassert>
#include <map>
#include <vector>

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/base/Algorithm.hpp>

namespace NetworKit {
    struct cc_result {
        node* components;
        node n_nodes;
        node* component_sizes;
        node n_components;
    };

/**
 * @ingroup components
 * Determines the connected components of an undirected graph.
 */
class ConnectedComponents : public Algorithm {
public:
    /**
     * Create ConnectedComponents class for Graph @a G.
     *
     * @param G The graph.
     */
    ConnectedComponents(const Graph& G);

    /**
     * This method determines the connected components for the graph given in the constructor.
     */
    void run() override;

    /**
     * Get the number of connected components.
     *
     * @return The number of connected components.
     */
    count numberOfComponents() const;

    /**
     * Get the the component in which node @a u is situated.
     *
     * @param[in]	u	The node whose component is asked for.
     */
    count componentOfNode(node u) const;


    /**
     * Get a Partition that represents the components.
     *
     * @return A partition representing the found components.
     */
    Partition getPartition() const;

    /**
     *Return the map from component to size
     */
    std::map<index, count> getComponentSizes() const;

    /**
     * @return Vector of components, each stored as (unordered) set of nodes.
     */
    std::vector<std::vector<node> > getComponents() const;

    /**
     * Constructs a new graph that contains only the nodes inside the largest
     * connected component.
     * @param G            The input graph.
     * @param compactGraph If true, the node ids of the output graph will be compacted
     * (i.e. re-numbered from 0 to n-1). If false, the node ids will not be changed.
     */
    static Graph extractLargestConnectedComponent(const Graph &G, bool compactGraph = false);
    static cc_result get_raw_partition(const Graph &G);

private:
    const Graph& G;
    Partition component;
    count numComponents;
};

inline count ConnectedComponents::componentOfNode(node u) const {
    assert (component[u] != none);
    assureFinished();
    return component[u];
}

inline count ConnectedComponents::numberOfComponents() const {
    assureFinished();
    return this->numComponents;
}

}


#endif // NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_HPP_
