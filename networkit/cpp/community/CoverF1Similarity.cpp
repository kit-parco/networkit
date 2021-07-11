// no-networkit-format

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/community/CoverF1Similarity.hpp>

namespace NetworKit {

CoverF1Similarity::CoverF1Similarity(const Graph &G, const Cover &C, const Cover &reference) : LocalCoverEvaluation(G, C), reference(&reference) {}

void CoverF1Similarity::run() {
    hasRun = false;
    Aux::SignalHandler handler;

    std::vector<std::vector<node>> Csets(C->upperBound());
    std::vector<count> referenceSizes(reference->upperBound(), 0);
    count numMemberships = 0;

    handler.assureRunning();

    // Calculate cluster sizes and explicitly store all clusters of C
    G->forNodes([&](node u) {
        for (index c : (*C)[u]) {
            Csets[c].push_back(u);
        }

        for (index c : (*reference)[u]) {
            ++referenceSizes[c];
        }

        numMemberships += (*C)[u].size();
    });

    handler.assureRunning();

    // Initialize result values
    unweightedAverage = 0;
    weightedAverage = 0;
    minimumValue = std::numeric_limits<double>::max();
    maximumValue = std::numeric_limits<double>::lowest();
    values.clear();
    values.resize(C->upperBound(), 0);

    count numClusters = 0;

    // The overlap of the currently considered cluster in C with each reference cluster
    std::vector<count> overlap(reference->upperBound(), 0);
    // The clusters of the reference that have a positive overlap
    std::vector<index> overlappingReference;

    for (index i = 0; i < C->upperBound(); ++i) {
        if (!Csets[i].empty()) {
            ++numClusters;

            for (node u : Csets[i]) {
                for (index s : (*reference)[u]) {
                    if (overlap[s] == 0) {
                        overlappingReference.push_back(s);
                    }

                    ++overlap[s];
                }
            }

            double bestF1 = 0;

            for (index s : overlappingReference) {
                count ol = overlap[s];
                assert(ol > 0);

                // Reset values
                overlap[s] = 0;

                double precision = ol * 1.0 / referenceSizes[s];
                double recall = ol * 1.0 / Csets[i].size();

                double f1 = 2 * (precision * recall) / (precision + recall);
                if (f1 > bestF1) {
                    bestF1 = f1;
                }
            }

            overlappingReference.clear();

            values[i] = bestF1;
            minimumValue = std::min(bestF1, minimumValue);
            maximumValue = std::max(bestF1, maximumValue);
            unweightedAverage += bestF1;
            weightedAverage += bestF1 * Csets[i].size();
        }
    }

    handler.assureRunning();

    unweightedAverage /= numClusters;
    weightedAverage /= numMemberships;

    hasRun = true;
}

}
