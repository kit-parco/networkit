/*
 * ConnectedComponents.hpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#ifndef NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_HPP_

#include <memory>

#include <networkit/components/ComponentDecomposition.hpp>

namespace NetworKit {

// pImpl
namespace ConnectedComponentsDetails {
template <bool>
class ConnectedComponentsImpl;
} // namespace ConnectedComponentsDetails

class ConnectedComponents final : public ComponentDecomposition {
public:
    /* Creates the ConnectedComponents class for graph @G.
     *
     * @param G The graph.
     */
    ConnectedComponents(const Graph &G);

    ~ConnectedComponents() override;

    /*
     * Computes the connected components of the input graph.
     */
    void run() override;

    /**
     * Constructs a new graph that contains only the nodes inside the largest
     * connected component.
     * @param G            The input graph.
     * @param compactGraph If true, the node ids of the output graph will be compacted
     * (i.e. re-numbered from 0 to n-1). If false, the node ids will not be changed.
     * @return The largest connected component of the input graph @a G.
     */
    static Graph extractLargestConnectedComponent(const Graph &G, bool compactGraph = false);

private:
    std::unique_ptr<ConnectedComponentsDetails::ConnectedComponentsImpl<false>> impl;
};

} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_HPP_
