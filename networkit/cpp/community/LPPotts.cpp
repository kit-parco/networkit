/*
 * LPPotts.cpp
 *
 * Created on: 2019-01-14
 * Author: Armin Wiebigke
 */

#include <omp.h>
#include <algorithm>
#include <cmath>

#include <networkit/community/LPPotts.hpp>
#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/UniformRandomSelector.hpp>
#include <networkit/auxiliary/SparseVector.hpp>
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/community/LPDegreeOrdered.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/auxiliary/Parallel.hpp>

namespace NetworKit {

LPPotts::LPPotts(const Graph &G, double alpha, count theta, count maxIterations,
                 bool parallelPropagation)
        : LPPotts(G, Partition(), alpha, theta, maxIterations, parallelPropagation) {
}

LPPotts::LPPotts(const Graph &G, const Partition &baseClustering, double alpha, count theta,
                 count maxIterations, bool parallel)
        : CommunityDetectionAlgorithm(G, baseClustering), alpha(alpha), updateThreshold(theta),
          maxIterations(maxIterations), parallel(parallel) {
}

void LPPotts::init() {
    numberOfThreads = parallel ? Aux::getMaxNumberOfThreads() : 1;
    index nodeIdBound = G->upperNodeIdBound();

    if (result.numberOfElements() == 0) {
        result = Partition(nodeIdBound);
        result.allToSingletons();
    }
    assert(result.numberOfElements() == nodeIdBound);

    if (updateThreshold == none) {
        updateThreshold = (count) (G->numberOfNodes() / 1e5);
        updateThreshold = std::max(updateThreshold, (count) 1);
    }
    globalLabelCounts.resize(nodeIdBound, 1);
    activeNodes.resize(nodeIdBound);
    G->forNodes([&](node u) { activeNodes[u] = true; });

    neighborLabelCountsPerThread.resize(numberOfThreads, SparseVector<count>(nodeIdBound, none));
    if (LPPotts::parallel) {
        globalLabelCountChangePerThread.resize(numberOfThreads,
                                               std::vector<int64_t>(nodeIdBound, 0));
        nextActiveNodesPerThread.resize(numberOfThreads, std::vector<uint8_t>(nodeIdBound, false));
    }
}

void LPPotts::run() {
    if (hasRun) {
        throw std::runtime_error("The algorithm has already run on the graph.");
    }
    init();
    runAlgorithm();
    hasRun = true;
}

void LPPotts::runAlgorithm() {
    const index nodeIdBound = G->upperNodeIdBound();
    count updatedNodesCount = G->numberOfNodes();
    iteration = 0;
    Partition secondPartition;
    std::vector<node> nodes;
    if (parallel) {
        secondPartition = result;
    } else {
        nodes = GraphTools::nodeSet(*G);
    }

    // Propagate labels
    while ((updatedNodesCount > updateThreshold) && (iteration < maxIterations)) {
        iteration += 1;
        updatedNodesCount = 0;
        DEBUG("[BEGIN] LabelPropagation: iteration #", iteration);

        if (parallel) {
            Aux::Timer timer{};
            timer.start();

#pragma omp parallel
            {
                count localUpdatedNodes = 0;
#pragma omp for nowait
                for (omp_index u = 0; u < static_cast<omp_index>(nodeIdBound); ++u) {
                    if (G->hasNode(u)) {
                        bool updated = evaluateNode(u, secondPartition);
                        if (updated)
                            ++localUpdatedNodes;
                    }
                }
#pragma omp atomic
                updatedNodesCount += localUpdatedNodes;
            }

            result = secondPartition;
            for (index tid = 0; tid < numberOfThreads; ++tid) {
#pragma omp parallel for schedule(static)
                for (omp_index i = 0; i < static_cast<omp_index>(nodeIdBound); ++i) {
                    // for label i
                    globalLabelCounts[i] += globalLabelCountChangePerThread[tid][i];
                    globalLabelCountChangePerThread[tid][i] = 0;
                    // for node i
                    activeNodes[i] |= nextActiveNodesPerThread[tid][i];
                    nextActiveNodesPerThread[tid][i] = false;
                }
            }

            timer.stop();
            timing.push_back(timer.elapsedMilliseconds());
            DEBUG("[DONE] LabelPropagation: iteration #", iteration, " - updated ",
                  updatedNodesCount, " labels, time spent: ", timer.elapsedTag());
        } else {
            std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());
            for (node u : nodes) {
                bool updated = evaluateNode(u, result);
                if (updated)
                    ++updatedNodesCount;
            }
        }
    } // end while
}

// Returns true if the label of the node was changed
bool LPPotts::evaluateNode(node u, Partition &nextPartition) {
    if (!activeNodes[u] || G->degree(u) == 0)
        return false;
    activeNodes[u] = false;

    label bestLabel = calculateBestLabel(u);

    label currentLabel = result.subsetOf(u);
    if (currentLabel != bestLabel) {
        updateLabel(u, nextPartition, currentLabel, bestLabel);
        return true;
    }
    return false;
}

void LPPotts::updateLabel(node u, Partition &nextPartition, label currentLabel, label newLabel) {
    index threadId = getThreadId();
    auto &globalLabelCountChange = parallel ? globalLabelCountChangePerThread[threadId]
                                            : globalLabelCounts;
    auto &nextActiveNodesThisThread = parallel ? nextActiveNodesPerThread[threadId]
                                               : activeNodes;
    nextPartition.moveToSubset(newLabel, u);
    --globalLabelCountChange[currentLabel];
    ++globalLabelCountChange[newLabel];
    nextActiveNodesThisThread[u] = true;
    G->forNeighborsOf(u, [&](node w) {
        nextActiveNodesThisThread[w] = true;
    });
}

label LPPotts::calculateBestLabel(node u) {
    // Count the labels of the neighbors
    SparseVector<count> &neighborLabelCounts = neighborLabelCountsPerThread[getThreadId()];
    assert(neighborLabelCounts.isClean());
    G->forNeighborsOf(u, [&](node w) {
        label neighborLabel = result.subsetOf(w);
        if (!neighborLabelCounts.indexIsUsed(neighborLabel)) {
            neighborLabelCounts.insert(neighborLabel, 0);
        }
        ++neighborLabelCounts[neighborLabel];
    });

    // Get best label
    Aux::UniformRandomSelector selector{};
    label bestLabel = none;
    double bestWeight = -std::numeric_limits<double>::max();
    for (label neighborLabel : neighborLabelCounts.insertedIndexes()) {
        count localCount = neighborLabelCounts[neighborLabel];
        double weight = localCount - alpha * (globalLabelCounts[neighborLabel] - localCount);
        if (weight > bestWeight) {
            bestWeight = weight;
            bestLabel = neighborLabel;
            selector.reset();
        } else if (weight == bestWeight) {
            // if multiple labels have the same weight, select one randomly
            if (selector.addElement())
                bestLabel = neighborLabel;
        }
    }
    assert(bestLabel != none);
    neighborLabelCounts.reset();
    return bestLabel;
}

index LPPotts::getThreadId() const {
    return parallel ? omp_get_thread_num() : 0;
}

std::string LPPotts::toString() const {
    std::stringstream strm;
    strm << "LPPotts";
    return strm.str();
}

void LPPotts::setUpdateThreshold(count th) {
    this->updateThreshold = th;
}

count LPPotts::numberOfIterations() {
    return this->iteration;
}

std::vector<count> LPPotts::getTiming() {
    return this->timing;
}

LPPottsFactory::LPPottsFactory(double alpha, count theta, count maxIterations,
                               bool parallelPropagation)
        : alpha(alpha), theta(theta), maxIterations(maxIterations),
          parallelPropagation(parallelPropagation) {
}

ClusteringFunction LPPottsFactory::getFunction() const {
    const double alphaCopy = alpha;
    const count thetaCopy = theta;
    const count maxIterationsCopy = maxIterations;
    const bool parallelPropagationCopy = parallelPropagation;
    return [alphaCopy, thetaCopy, maxIterationsCopy, parallelPropagationCopy](const Graph &G) {
        LPPotts lpPotts(G, alphaCopy, thetaCopy, maxIterationsCopy, parallelPropagationCopy);
        lpPotts.run();
        return lpPotts.getPartition();
    };
}
} /* namespace NetworKit */
