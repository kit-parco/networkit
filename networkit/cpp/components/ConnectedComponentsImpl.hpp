
#ifndef NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_GENERAL_HPP_
#define NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_GENERAL_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {
namespace ConnectedComponentsDetails {

/**
 * @ingroup components
 * Determines the (weakly) connected components of a graph.
 */
template <bool WeaklyCC = false>
class ConnectedComponentsImpl final : public Algorithm {
public:
    /* Create the ConnectedComponentsImpl class for graph @G.
     *
     * @param G The graph.
     */
    ConnectedComponentsImpl(const Graph &G, Partition &components);

    /*
     * Compute the (weakly) connected components of the input graph.
     */
    void run() override;

    /**
     * Constructs a new graph that contains only the nodes inside the largest
     * (weakly) connected component.
     * @param G            The input graph.
     * @param compactGraph If true, the node ids of the output graph will be compacted
     * (i.e. re-numbered from 0 to n-1). If false, the node ids will not be changed.
     * @return The largest (weakly) connected component of the input graph @a G.
     */
    static Graph extractLargestConnectedComponent(const Graph &G, bool compactGraph);

private:
    const Graph *G;
    Partition *componentPtr;
};

} // namespace ConnectedComponentsDetails
} // namespace NetworKit

#endif /* ifndef NETWORKIT_COMPONENTS_CONNECTED_COMPONENTS_GENERAL_HPP_ */
