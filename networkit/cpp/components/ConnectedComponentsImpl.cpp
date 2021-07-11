
#include "ConnectedComponentsImpl.hpp"

#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {
namespace ConnectedComponentsDetails {

template <bool WeaklyCC>
ConnectedComponentsImpl<WeaklyCC>::ConnectedComponentsImpl(const Graph &G, Partition &components)
    : G(&G), componentPtr(&components) {

    if (!WeaklyCC && G.isDirected())
        throw std::runtime_error(
            "Error, connected components of directed graphs cannot be "
            "computed, use StronglyConnectedComponents or WeaklyConnectedComponents instead.");

    if (WeaklyCC && !G.isDirected())
        throw std::runtime_error("Error, weakly connected components of an undirected graph cannot "
                                 "be computed, use ConnectedComponents instead.");
}

template <bool WeaklyCC>
void ConnectedComponentsImpl<WeaklyCC>::run() {
    index nComponents = 0;
    count visitedNodes = 0;

    std::queue<node> q;
    auto &component = *componentPtr;
    component.reset(G->upperNodeIdBound(), none);

    // perform breadth-first searches
    for (const node u : G->nodeRange()) {
        if (component[u] != none)
            continue;

        component.setUpperBound(nComponents + 1);
        const index c = nComponents;

        q.push(u);
        component[u] = c;

        do {
            const node v = q.front();
            q.pop();
            ++visitedNodes;

            const auto visitNeighbor = [&](node w) -> void {
                if (component[w] == none) {
                    q.push(w);
                    component[w] = c;
                }
            };

            // enqueue neighbors, set component
            G->forNeighborsOf(v, visitNeighbor);
            if (WeaklyCC)
                G->forInNeighborsOf(v, visitNeighbor);

        } while (!q.empty());

        ++nComponents;

        if (visitedNodes == G->numberOfNodes())
            break;
    }

    hasRun = true;
}

template <bool WeaklyCC>
Graph ConnectedComponentsImpl<WeaklyCC>::extractLargestConnectedComponent(const Graph &G,
                                                                          bool compactGraph) {
    if (G.isEmpty())
        return G;

    Partition component;
    ConnectedComponentsImpl<WeaklyCC> cc(G, component);
    cc.run();
    const auto compSizes = component.subsetSizeMap();
    if (compSizes.size() == 1) {
        if (compactGraph)
            return GraphTools::getCompactedGraph(G, GraphTools::getContinuousNodeIds(G));
        return G;
    }

    const auto largestCCIndex =
        std::max_element(compSizes.begin(), compSizes.end(),
                         [](const std::pair<index, count> &x, const std::pair<index, count> &y) {
                             return x.second < y.second;
                         })
            ->first;

    const auto nodesInLCC = component.getMembers(largestCCIndex);
    return GraphTools::subgraphFromNodes(G, nodesInLCC.begin(), nodesInLCC.end(), compactGraph);
}

template class ConnectedComponentsImpl<false>;
template class ConnectedComponentsImpl<true>;

} // namespace ConnectedComponentsDetails
} // namespace NetworKit
