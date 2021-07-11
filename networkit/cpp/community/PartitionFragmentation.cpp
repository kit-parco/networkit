// no-networkit-format
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/community/PartitionFragmentation.hpp>
#include <networkit/community/PartitionIntersection.hpp>
#include <networkit/components/ConnectedComponents.hpp>

void NetworKit::PartitionFragmentation::run() {
    hasRun = false;

    Aux::SignalHandler handler;

    ConnectedComponents cc(*G);
    cc.run();

    handler.assureRunning();

    Partition ccP = cc.getPartition();

    handler.assureRunning();

    Partition ints = PartitionIntersection().calculate(*P, ccP);

    handler.assureRunning();

    std::vector<count> intsSizes(ints.upperBound());
    std::vector<count> PSizes(P->upperBound());

    handler.assureRunning();

    G->forNodes([&](node u) {
        ++intsSizes[ints[u]];
        ++PSizes[(*P)[u]];
    });

    handler.assureRunning();

    values.clear();
    values.resize(P->upperBound(), std::numeric_limits< double >::max());

    handler.assureRunning();

    G->forNodes([&](node u) {
        values[(*P)[u]] = std::min(values[(*P)[u]], 1.0 - intsSizes[ints[u]] * 1.0 / PSizes[(*P)[u]]);
    });

    handler.assureRunning();

    maximumValue = 0.0;
    minimumValue = std::numeric_limits<double>::max();
    unweightedAverage = 0.0;
    weightedAverage = 0.0;

    count numSubsets = 0;

    for (index i = 0; i < P->upperBound(); ++i) {
        if (values[i] < std::numeric_limits< double >::max()) {
            maximumValue = std::max(values[i], maximumValue);
            minimumValue = std::min(values[i], minimumValue);
            unweightedAverage += values[i];
            weightedAverage += values[i] * PSizes[i];

            ++numSubsets;
        }
    }

    handler.assureRunning();

    unweightedAverage /= numSubsets;
    weightedAverage /= G->numberOfNodes();

    hasRun = true;
}

