// no-networkit-format
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/community/IsolatedInterpartitionExpansion.hpp>

void NetworKit::IsolatedInterpartitionExpansion::run() {
    hasRun = false;

    Aux::SignalHandler handler;

    values.clear();
    values.resize(P->upperBound(), 0);

    handler.assureRunning();

    G->forEdges([&](node u, node v, edgeweight w) {
        if ((*P)[u] != (*P)[v]) {
            values[(*P)[u]] += w;
            values[(*P)[v]] += w;
        }
    });

    handler.assureRunning();

    std::vector<count> clusterSize(P->upperBound(), 0);
    G->forNodes([&](node u) {
        ++clusterSize[(*P)[u]];
    });

    handler.assureRunning();

    weightedAverage = 0;
    unweightedAverage = 0;
    maximumValue = std::numeric_limits<double>::lowest();
    minimumValue = std::numeric_limits<double>::max();

    count n = G->numberOfNodes();
    count c = 0;

    for (index i = 0; i < P->upperBound(); ++i) {
        if (clusterSize[i] > 0) {
            edgeweight cond = values[i] / std::min(clusterSize[i], n - clusterSize[i]);
            unweightedAverage += cond;
            weightedAverage += clusterSize[i] * cond;
            maximumValue = std::max(maximumValue, cond);
            minimumValue = std::min(minimumValue, cond);
            ++c;
        }
    }

    handler.assureRunning();

    unweightedAverage /= c;
    weightedAverage /= G->numberOfNodes();

    hasRun = true;
}
