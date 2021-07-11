// no-networkit-format
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/community/IntrapartitionDensity.hpp>

void NetworKit::IntrapartitionDensity::run() {
    hasRun = false;

    Aux::SignalHandler handler;

    minimumValue = std::numeric_limits< double >::max();
    maximumValue  = std::numeric_limits< double >::lowest();
    unweightedAverage = 0;
    weightedAverage = 0;
    values.clear();

    std::vector<count> clusterSizes(P->upperBound(), 0);
    std::vector<count> intraEdges(P->upperBound(), 0);

    handler.assureRunning();

    G->forEdges([&](node u, node v) {
        if ((*P)[u] == (*P)[v]) {
            ++intraEdges[(*P)[u]];
        }
    });

    handler.assureRunning();

    G->forNodes([&](node u) {
        ++clusterSizes[(*P)[u]];
    });

    handler.assureRunning();

    count numClusters = 0;
    count intraEdgesSum = 0;
    count possibleIntraEdgesSum = 0;

    values.resize(P->upperBound(), 0);

    for (index i = 0; i < clusterSizes.size(); ++i) {
        if (clusterSizes[i] > 0) {
            double id = 1;
            count possibleEdges = clusterSizes[i] * (clusterSizes[i]-1) / 2;
            if (possibleEdges > 0) {
                id = intraEdges[i] * 1.0 / possibleEdges;
            }

            values[i] = id;

            unweightedAverage += id;
            weightedAverage += id * clusterSizes[i];
            minimumValue = std::min(id, minimumValue);
            maximumValue = std::max(id, maximumValue);

            possibleIntraEdgesSum += possibleEdges;
            intraEdgesSum += intraEdges[i];
            ++numClusters;
        } else {
            clusterSizes[i] = none;
        }
    }

    handler.assureRunning();

    unweightedAverage /= numClusters;
    weightedAverage /= G->numberOfNodes();

    globalValue = intraEdgesSum * 1.0 / possibleIntraEdgesSum;
    hasRun = true;
}
