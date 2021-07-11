
#include <networkit/components/ComponentDecomposition.hpp>

namespace NetworKit {

ComponentDecomposition::ComponentDecomposition(const Graph &G)
    : G(&G), component(G.upperNodeIdBound()) {
    hasRun = false;
}

count ComponentDecomposition::numberOfComponents() const {
    assureFinished();
    return component.upperBound();
}

count ComponentDecomposition::componentOfNode(node u) const {
    assureFinished();
    assert(component[u] != none);
    return component[u];
}

const Partition &ComponentDecomposition::getPartition() const {
    assureFinished();
    return component;
}

std::vector<std::vector<node>> ComponentDecomposition::getComponents() const {
    assureFinished();
    std::vector<std::vector<node>> result(numberOfComponents());
    G->forNodes([&](node u) { result[component[u]].push_back(u); });
    return result;
}

std::map<index, count> ComponentDecomposition::getComponentSizes() const {
    assureFinished();
    return component.subsetSizeMap();
}

} // namespace NetworKit
