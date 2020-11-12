#include <vector>

#include "DynApproxBetweennessImpl.hpp"
#include <networkit/auxiliary/Log.hpp>
#include <networkit/centrality/DynApproxBetweenness.hpp>
#include <networkit/dynamics/GraphEvent.hpp>

// networkit-format
namespace NetworKit {

DynApproxBetweenness::DynApproxBetweenness(const Graph &G, double epsilon, double delta,
                                           bool storePredecessors, double universalConstant)
    : Centrality(G, true),
      impl(new DynApproxBetweennessImpl(G, epsilon, delta, storePredecessors, universalConstant)) {}

DynApproxBetweenness::~DynApproxBetweenness() = default;

void DynApproxBetweenness::run() {
    impl->run();
    hasRun = true;
}

void DynApproxBetweenness::update(GraphEvent e) {
    impl->update(e);
}

void DynApproxBetweenness::updateBatch(const std::vector<GraphEvent> &e) {
    impl->updateBatch(e);
}

count DynApproxBetweenness::getNumberOfSamples() {
    return impl->getNumberOfSamples();
}
} // namespace NetworKit
