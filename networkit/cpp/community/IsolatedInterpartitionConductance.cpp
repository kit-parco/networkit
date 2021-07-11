// no-networkit-format
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/community/IsolatedInterpartitionConductance.hpp>

void NetworKit::IsolatedInterpartitionConductance::run() {
    hasRun = false;

    Aux::SignalHandler handler;

    values.clear();
    values.resize(P->upperBound(), 0);

    handler.assureRunning();

    std::vector<edgeweight> clusterVolume(P->upperBound(), 0);
    edgeweight totalVolume = 0;
    G->forEdges([&](node u, node v, edgeweight w) {
        if ((*P)[u] != (*P)[v]) {
            values[(*P)[u]] += w;
            values[(*P)[v]] += w;
        }

        clusterVolume[(*P)[u]] += w;
        clusterVolume[(*P)[v]] += w;
        totalVolume += 2 * w;
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

    count c = 0;

    for (index i = 0; i < P->upperBound(); ++i) {
        if (clusterSize[i] > 0) {
            edgeweight cond = 0;
            auto denominator = std::min(clusterVolume[i], totalVolume - clusterVolume[i]);

            if (denominator > 0) {
                cond = values[i] / denominator;
            }

            values[i] = cond;
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
