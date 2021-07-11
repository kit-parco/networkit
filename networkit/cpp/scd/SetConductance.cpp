#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/scd/SetConductance.hpp>

namespace NetworKit {
SetConductance::SetConductance(const Graph &g, const std::set<node> &community)
    : g(&g), community(&community) {
    if (g.isDirected()) {
        throw std::runtime_error("SetConductance only supports undirected graphs.");
    }
}

void SetConductance::run() {
    hasRun = false;

    Aux::SignalHandler handler;

    double cut = 0;
    double comVolume = 0;

    for (node u : *community) {
        if (g->hasNode(u)) {
            g->forEdgesOf(u, [&](node, node v, edgeweight ew) {
                if (community->count(v) == 0) {
                    cut += ew;
                }
                comVolume += ew;

                // count loops twice
                if (u == v) {
                    comVolume += ew;
                }
            });
        }
    }

    double totalVolume = g->totalEdgeWeight() * 2;

    double restVolume = totalVolume - comVolume;

    conductance = 1.0;

    if (comVolume > 0 && restVolume > 0) {
        conductance = cut / std::min(comVolume, restVolume);
    }

    hasRun = true;
}

double SetConductance::getConductance() const {
    assureFinished();
    return conductance;
}

} // namespace NetworKit
