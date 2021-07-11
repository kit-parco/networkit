// no-networkit-format
#include <atomic>
#include <memory>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/community/PartitionHubDominance.hpp>

void NetworKit::PartitionHubDominance::run() {
    hasRun = false;

    Aux::SignalHandler handler;

    std::unique_ptr<std::atomic<count>[]> maxInternalDeg(new std::atomic<count>[P->upperBound()]{});
    std::vector<count> clusterSizes(P->upperBound(), 0);

    handler.assureRunning();

    G->balancedParallelForNodes([&](node u) {
        index c = (*P)[u];

        if (c != none) {
            count internalDeg = 0;
            G->forNeighborsOf(u, [&](node v) {
                if ((*P)[v] == c) {
                    internalDeg++;
                }
            });

            Aux::Parallel::atomic_max(maxInternalDeg[c], internalDeg);

            #pragma omp atomic
            ++clusterSizes[c];
        }
    });

    handler.assureRunning();

    count numClusters = 0;
    weightedAverage = 0;
    unweightedAverage = 0;
    maximumValue = std::numeric_limits<double>::lowest();
    minimumValue = std::numeric_limits<double>::max();
    values.clear();
    values.resize(P->upperBound(), 0);

    for (index i = 0; i < P->upperBound(); ++i) {
        if (clusterSizes[i] > 0) {
            ++numClusters;

            double dominance = 1;
            if (clusterSizes[i] > 1) {
                dominance = static_cast<double>(maxInternalDeg[i]) * 1.0 / static_cast<double>(clusterSizes[i] - 1);
            }

            values[i] = dominance;
            unweightedAverage += dominance;
            weightedAverage = dominance * clusterSizes[i];

            maximumValue = std::max(dominance, maximumValue);
            minimumValue = std::min(dominance, minimumValue);
        }
    }

    handler.assureRunning();

    unweightedAverage /= numClusters;
    weightedAverage /= G->numberOfNodes();
    hasRun = true;
}
