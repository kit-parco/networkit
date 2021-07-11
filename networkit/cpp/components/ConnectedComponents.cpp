/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#include "ConnectedComponentsImpl.hpp"

#include <networkit/components/ConnectedComponents.hpp>

namespace NetworKit {

ConnectedComponents::ConnectedComponents(const Graph &G)
    : ComponentDecomposition(G),
      impl(new ConnectedComponentsDetails::ConnectedComponentsImpl<false>{G, component}) {}

ConnectedComponents::~ConnectedComponents() = default;

void ConnectedComponents::run() {
    impl->run();
    hasRun = true;
}

Graph ConnectedComponents::extractLargestConnectedComponent(const Graph &G, bool compactGraph) {
    return ConnectedComponentsDetails::ConnectedComponentsImpl<
        false>::extractLargestConnectedComponent(G, compactGraph);
}

} // namespace NetworKit
