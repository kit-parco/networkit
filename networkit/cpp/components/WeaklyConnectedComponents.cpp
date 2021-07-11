/*
 * WeaklyConnectedComponents.cpp
 *
 *  Created on: June 20, 2017
 *      Author: Eugenio Angriman
 */

#include "ConnectedComponentsImpl.hpp"

#include <networkit/components/WeaklyConnectedComponents.hpp>

namespace NetworKit {

WeaklyConnectedComponents::WeaklyConnectedComponents(const Graph &G)
    : ComponentDecomposition(G),
      impl(new ConnectedComponentsDetails::ConnectedComponentsImpl<true>{G, component}) {}

WeaklyConnectedComponents::~WeaklyConnectedComponents() = default;

void WeaklyConnectedComponents::run() {
    impl->run();
    hasRun = true;
}

Graph WeaklyConnectedComponents::extractLargestWeaklyConnectedComponent(const Graph &G,
                                                                        bool compactGraph) {
    return ConnectedComponentsDetails::ConnectedComponentsImpl<
        true>::extractLargestConnectedComponent(G, compactGraph);
}

} // namespace NetworKit
