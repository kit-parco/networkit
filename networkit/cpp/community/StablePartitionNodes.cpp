// no-networkit-format
#include <limits>
#include <map>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/community/StablePartitionNodes.hpp>

void NetworKit::StablePartitionNodes::run() {
    hasRun = false;

    Aux::SignalHandler handler;

    stableMarker.clear();
    stableMarker.resize(G->upperNodeIdBound(), true);
    values.clear();

    handler.assureRunning();

    // first determine which nodes are stable
    G->balancedParallelForNodes([&](node u) {
        if (G->degree(u) > 0) { // we consider isolated nodes to be stable.
            std::map<index, count> labelWeights;
            G->forNeighborsOf(u, [&](node v, edgeweight ew) {
                labelWeights[(*P)[v]] += ew;
            });

            index ownLabel = (*P)[u];
            double ownWeight = labelWeights[ownLabel];

            if (ownWeight == 0) {
                stableMarker[u] = false;
            } else {
                for (auto lw : labelWeights) {
                    if (lw.first != ownLabel && lw.second >= ownWeight) {
                        stableMarker[u] = false;
                        break;
                    }
                }
            }
        }
    });

    handler.assureRunning();

    values.resize(P->upperBound(), 0);
    std::vector<count> partitionSizes(P->upperBound(), 0);
    count stableCount = 0;

    // collect how many nodes are stable in which partition
    G->forNodes([&](node u) {
        ++partitionSizes[(*P)[u]];
        values[(*P)[u]] += stableMarker[u];
        stableCount += stableMarker[u];
    });

    count numClusters = 0;
    unweightedAverage = 0;
    minimumValue = std::numeric_limits<double>::max();
    maximumValue = std::numeric_limits<double>::lowest();

    // calculate all average/max/min-values
    for (index i = 0; i < P->upperBound(); ++i) {
        if (partitionSizes[i] > 0) {
            values[i] /= partitionSizes[i];
            unweightedAverage += values[i];
            minimumValue = std::min(minimumValue, values[i]);
            maximumValue = std::max(maximumValue, values[i]);
            ++numClusters;
        }
    }

    unweightedAverage /= numClusters;
    weightedAverage = stableCount * 1.0 / G->numberOfNodes();

    handler.assureRunning(); // make sure we do not ignore the signal sent by the user

    hasRun = true;
}
