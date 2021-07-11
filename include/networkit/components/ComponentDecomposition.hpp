/*
 * ComponentDecomposition.hpp
 *
 *  Created on: 17.12.2020
 *      Author: cls,
 *              Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_COMPONENTS_COMPONENT_DECOMPOSITION_HPP_
#define NETWORKIT_COMPONENTS_COMPONENT_DECOMPOSITION_HPP_

#include <map>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup components
 * Abstract class for algorithms that compute the components of a graph.
 */
class ComponentDecomposition : public Algorithm {
public:
    /**
     * Constructs the ComponentDecomposition class for the given Graph @a G.
     *
     * @param G The graph.
     */
    ComponentDecomposition(const Graph &G);

    /**
     * Get the number of connected components.
     *
     * @return The number of connected components.
     */
    count numberOfComponents() const;

    /**
     * Get the component in which node @a u is situated.
     *
     * @param[in] u The node whose component is asked for.
     */
    count componentOfNode(node u) const;

    /**
     * Get a Partition that represents the components.
     *
     * @return A partition representing the found components.
     */
    const Partition &getPartition() const;

    /**
     * Return the map from component to size
     */
    std::map<index, count> getComponentSizes() const;

    /**
     * @return Vector of components, each stored as (unordered) set of nodes.
     */
    std::vector<std::vector<node>> getComponents() const;

protected:
    const Graph *G;
    Partition component;
};

} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_COMPONENT_DECOMPOSITION_HPP_
