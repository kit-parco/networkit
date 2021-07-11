/*
 * WeaklyConnectedComponents.hpp
 *
 *  Created on: June 20, 2017
 *      Author: Eugenio Angriman
 */

#ifndef NETWORKIT_COMPONENTS_WEAKLY_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_WEAKLY_CONNECTED_COMPONENTS_HPP_

#include <memory>

#include <networkit/components/ComponentDecomposition.hpp>

namespace NetworKit {

// pImpl
namespace ConnectedComponentsDetails {
template <bool>
class ConnectedComponentsImpl;
} // namespace ConnectedComponentsDetails

/**
 * @ingroup components
 * Determines the weakly connected components of a directed graph.
 */
class WeaklyConnectedComponents final : public ComponentDecomposition {
public:
    /**
     * Create WeaklyConnectedComponents class for Graph @a G.
     *
     * @param G The graph.
     */
    WeaklyConnectedComponents(const Graph &G);

    ~WeaklyConnectedComponents() override;

    /*
     * Computes the weakly connected components of the input graph.
     */
    void run() override;

    /**
     * Constructs a new graph that contains only the nodes inside the largest
     * weakly connected component.
     * @param G            The input graph.
     * @param compactGraph If true, the node ids of the output graph will be compacted
     * (i.e. re-numbered from 0 to n-1). If false, the node ids will not be changed.
     * @return The largest weakly connected component of the input graph @a G.
     */
    static Graph extractLargestWeaklyConnectedComponent(const Graph &G, bool compactGraph = false);

private:
    std::unique_ptr<ConnectedComponentsDetails::ConnectedComponentsImpl<true>> impl;
};
} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_WEAKLY_CONNECTED_COMPONENTS_HPP_
