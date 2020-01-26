/*
 * MergeCommunities.cpp
 *
 * Created: 2019-09-26
 * Author: Armin Wiebigke
 */

#include <unordered_map>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/community/cleanup/MergeCommunities.hpp>
#include <networkit/graph/GraphBuilder.hpp>

namespace NetworKit {

MergeCommunities::MergeCommunities(const Graph &graph,
                                   std::vector<std::vector<node>> discardedCommunities,
                                   StochasticDistributionCalculator &stochasticDistribution,
                                   double significanceThreshold,
                                   double scoreThreshold,
                                   double minOverlapRatio,
                                   count maxCommunitySize)
        : graph(graph),
          discardedCommunities(std::move(discardedCommunities)),
          stochasticDistribution(stochasticDistribution),
          significanceCalculator({stochasticDistribution}),
          significanceThreshold(significanceThreshold),
          scoreThreshold(scoreThreshold),
          minOverlapRatio(minOverlapRatio),
          maxCommunitySize(maxCommunitySize) {
}

void MergeCommunities::run() {
    Aux::Timer timer{};
    timer.start();
    auto printTime = [&](const std::string &name) {
        timer.stop();
        INFO(name, " took ", (double) timer.elapsedMilliseconds() / 1000, "s");
        timer.start();
    };
    count iterations = 2;
    for (count i = 0; i < iterations; ++i) {
        createDiscardedCommunitiesGraph();
        printTime("Creating the discarded communities graph");
        tryToMergeCommunities();
        printTime("Local move");
        checkMergedCommunities();
        printTime("Clean merged communities");
    }
    hasRun = true;
}

void MergeCommunities::createDiscardedCommunitiesGraph() {
    count numDiscardedCommunities = discardedCommunities.size();

    std::vector<std::vector<index>> nodeMemberships(graph.upperNodeIdBound());
    for (index commId = 0; commId < discardedCommunities.size(); ++commId) {
        for (node u : discardedCommunities[commId]) {
            nodeMemberships[u].push_back(commId);
        }
    }

    outgoingGroupStubs.clear();
    outgoingGroupStubs.resize(numDiscardedCommunities);
    totalGroupStubs.clear();
    totalGroupStubs.resize(numDiscardedCommunities);
    totalStubs = 0;
    GraphBuilder builder(numDiscardedCommunities, true, false);
#pragma omp parallel
    {
        count localStubs = 0;
        SparseVector<count> neighborCommunities(graph.upperNodeIdBound(), 0);

        auto numComms = static_cast<omp_index>(discardedCommunities.size());
#pragma omp for schedule(dynamic, 10) nowait
        for (omp_index commId = 0; commId < numComms; ++commId) {
            for (node u : discardedCommunities[commId]) {
                graph.forNeighborsOf(u, [&](node v) {
                    for (index comm2 : nodeMemberships[v]) {
                        if (!neighborCommunities.indexIsUsed(comm2)) {
                            neighborCommunities.insert(comm2, 1);
                        } else {
                            neighborCommunities[comm2] += 1;
                        }
                    }
                });
            }

            for (index neighborCommunity : neighborCommunities.insertedIndexes()) {
                count weight = neighborCommunities[neighborCommunity];
                totalGroupStubs[commId] += weight;
                if (commId != neighborCommunity) {
                    outgoingGroupStubs[commId] += weight;
                } else {
                    // count loops twice
                    totalGroupStubs[commId] += weight;
                }
                builder.addHalfEdge(commId, neighborCommunity, weight);
            }

            localStubs += totalGroupStubs[commId];
            neighborCommunities.reset();
        }

#pragma omp atomic
        totalStubs += localStubs;
    }

    discardedCommunitiesGraph = builder.toGraph(false);

    INFO("Total stubs of discarded graph: ", totalStubs, ", number of edges in discarded graph: ",
         discardedCommunitiesGraph.numberOfEdges(), ", number of edges in original: ",
         graph.numberOfEdges());
    stochasticDistribution.increaseMaxValueTo(totalStubs + numDiscardedCommunities);
}

void MergeCommunities::tryToMergeCommunities() {
    mergedCommunities = Partition(discardedCommunitiesGraph.upperNodeIdBound());
    mergedCommunities.allToSingletons();
    SparseVector<edgeweight> neighborWeights(discardedCommunitiesGraph.upperNodeIdBound(), 0.0);

    count maxIterations = 20;
    for (count i = 0; i < maxIterations; ++i) {
        DEBUG("Merge communities iteration ", i);
        count nodesChanged = 0;
        discardedCommunitiesGraph.forNodesInRandomOrder([&](node u) {
            bool wasMoved = tryLocalMove(u, neighborWeights);
            if (wasMoved)
                ++nodesChanged;
        });
        if (nodesChanged == 0)
            break;
    }
    INFO("Merged ", mergedCommunities.numberOfElements(), " discarded communities into ",
         mergedCommunities.numberOfSubsets(), " communities");
}

bool MergeCommunities::tryLocalMove(node u, SparseVector<edgeweight> &neighborWeights) {
    edgeweight loopCount = 0;
    edgeweight degree = 0;
    discardedCommunitiesGraph.forEdgesOf(u, [&](node, node v, edgeweight weight) {
        auto neighborCommunity = mergedCommunities.subsetOf(v);
        degree += weight;
        if (v == u) { // count loops twice
            degree += weight;
            loopCount += weight;
        }
        if (!neighborWeights.indexIsUsed(neighborCommunity)) {
            neighborWeights.insert(neighborCommunity, weight);
        } else {
            neighborWeights[neighborCommunity] += weight;
        }
    });

    index currentCommunity = mergedCommunities.subsetOf(u);
    index bestNeighborCommunity = currentCommunity;
    double bestScore = 1.;
    count stubsIntoCurrent = 0;
    count stubsIntoNew = 0;
    for (index neighborCommunity : neighborWeights.insertedIndexes()) {
        double numEdgesIntoCommunity = neighborWeights[neighborCommunity];
        double score = 1.;
        count externalStubs = totalStubs - totalGroupStubs[neighborCommunity];
        if (neighborCommunity == currentCommunity) {
            // Calculate r-Score as if the node was not in the community
            numEdgesIntoCommunity -= loopCount;
            stubsIntoCurrent = numEdgesIntoCommunity;
            count outgoingFromNode = degree - 2 * loopCount - numEdgesIntoCommunity;
            count outgoingStubs =
                    outgoingGroupStubs[currentCommunity] - outgoingFromNode + numEdgesIntoCommunity;
            externalStubs += degree;
            bool communityIsEmpty = outgoingStubs != 0;
            if (!communityIsEmpty) {
                score = significanceCalculator.rScore(degree, numEdgesIntoCommunity,
                                                      outgoingStubs, externalStubs);
            }
        } else {
            score = significanceCalculator.rScore(degree, numEdgesIntoCommunity,
                                                  outgoingGroupStubs[neighborCommunity],
                                                  externalStubs);
        }
        if (score < bestScore) {
            stubsIntoNew = numEdgesIntoCommunity;
            bestScore = score;
            bestNeighborCommunity = neighborCommunity;
        }
    }

    neighborWeights.reset();

    if (bestNeighborCommunity != currentCommunity) {
        // Move node into new community
        mergedCommunities.moveToSubset(bestNeighborCommunity, u);
        count degreeWithoutLoop = degree - 2 * loopCount;
        outgoingGroupStubs[currentCommunity] += -degreeWithoutLoop + 2 * stubsIntoCurrent;
        outgoingGroupStubs[bestNeighborCommunity] += degreeWithoutLoop - 2 * stubsIntoNew;
        totalGroupStubs[currentCommunity] -= degree;
        totalGroupStubs[bestNeighborCommunity] += degree;
        return true;
    }
    return false;
}

void MergeCommunities::checkMergedCommunities() {
    // Store number of communities as this is not an O(1) lookup
    const count numMergedCommunities = mergedCommunities.numberOfSubsets();
    INFO("Check significance of ", numMergedCommunities, " merged communities");
    count skippedCommunities = 0;

    std::vector<index> partitionMap(mergedCommunities.upperBound(), none);
    std::vector<std::vector<node>> mergedCommunitiesSubsets;
    mergedCommunitiesSubsets.reserve(numMergedCommunities);
    mergedCommunities.forEntries([&](index e, index s) {
        if (partitionMap[s] == none) {
            partitionMap[s] = mergedCommunitiesSubsets.size();
            mergedCommunitiesSubsets.emplace_back();
        }
        mergedCommunitiesSubsets[partitionMap[s]].push_back(e);
    });

#pragma omp parallel
    {
        SingleCommunityCleanUp singleCommunityCleanUp(graph, stochasticDistribution, scoreThreshold,
                                                      significanceThreshold, minOverlapRatio);
        SparseVector<bool> mergedCommunity(graph.upperNodeIdBound());
#pragma omp for schedule(dynamic, 1)
        for (omp_index i = 0; i < static_cast<omp_index>(mergedCommunitiesSubsets.size()); ++i) {
            const std::vector<node> &communitiesToMerge = mergedCommunitiesSubsets[i];
            DEBUG("Clean merged community ", i, "/", numMergedCommunities);
            if (communitiesToMerge.size() == 1)
                continue;

            for (index communityId : communitiesToMerge) {
                const auto &community = discardedCommunities[communityId];
                if (mergedCommunity.size() + community.size() > maxCommunitySize) {
                    discardedCommunities[communityId].clear();
                    continue;
                }
                for (node u : community) {
                    if (!mergedCommunity.indexIsUsed(u)) {
                        mergedCommunity.insert(u, true);
                    }
                }
                discardedCommunities[communityId].clear();
            }
            if (mergedCommunity.size() > maxCommunitySize) {
#pragma omp atomic
                ++skippedCommunities;
            } else {
                assert(stochasticDistribution.maxValue() >=
                       2 * graph.totalEdgeWeight() + graph.numberOfNodes());
                const std::vector<index>& nodesOfMergedCommunity = mergedCommunity.insertedIndexes();
                std::vector<node> cleanedCommunity = singleCommunityCleanUp.clean(
                        nodesOfMergedCommunity);
                if (!cleanedCommunity.empty()) {
#pragma omp critical
                    cleanedCommunities.emplace_back(std::move(cleanedCommunity));
                } else {
                    discardedCommunities[communitiesToMerge.front()] = nodesOfMergedCommunity;
                }
            }

            mergedCommunity.reset();
        }
    }

    auto new_end = std::remove_if(discardedCommunities.begin(), discardedCommunities.end(),
                                  [](const std::vector<node> &c) { return c.empty(); });
    discardedCommunities.erase(new_end, discardedCommunities.end());
    INFO("Skipped ", skippedCommunities, " large communities (max size ", maxCommunitySize, ")");
}

const std::vector<std::vector<node>> &MergeCommunities::getCleanedCommunities() {
    return cleanedCommunities;
}

std::string MergeCommunities::toString() const {
    return "MergeCommunities";
}

bool MergeCommunities::isParallel() const {
    return true;
}

} /* namespace NetworKit */
