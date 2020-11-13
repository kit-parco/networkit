#include <vector>

#include "DynApproxBetweennessImpl.hpp"
#include "DynApproxBetweennessImplDir.hpp"
#include <networkit/auxiliary/Log.hpp>
#include <networkit/centrality/DynApproxBetweenness.hpp>
#include <networkit/dynamics/GraphEvent.hpp>

// networkit-format
namespace NetworKit {

DynApproxBetweenness::DynApproxBetweenness(const Graph &G, double epsilon, double delta,
                                           bool storePredecessors, double universalConstant)
    : Centrality(G, true),
      impl(new DynApproxBetweennessImpl(G, epsilon, delta, storePredecessors, universalConstant)) {}

DynApproxBetweenness::DynApproxBetweenness(const Graph &G, double epsilon, double delta,
                                           bool storePredecessors)
    : Centrality(G, true),
      implDir(new DynApproxBetweennessImplDir(G, epsilon, delta, storePredecessors)) {}

DynApproxBetweenness::~DynApproxBetweenness() = default;

void DynApproxBetweenness::run() {
    if (G.isDirected()) {
        implDir->run();
    } else {
        impl->run();
    }
    hasRun = true;
}

void DynApproxBetweenness::update(GraphEvent e) {
    if (G.isDirected()) {
        implDir->update(e);
    } else {
        impl->update(e);
    }
}

void DynApproxBetweenness::updateBatch(const std::vector<GraphEvent> &e) {
    if (G.isDirected()) {
        implDir->updateBatch(e);
    } else {
        impl->updateBatch(e);
    }
}

count DynApproxBetweenness::getNumberOfSamples() {
    if (G.isDirected()) {
        return implDir->getNumberOfSamples();
    } else {
        return impl->getNumberOfSamples();
    }
}
} // namespace NetworKit
