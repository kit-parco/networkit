// no-networkit-format
#include <algorithm>
#include <numeric>
#include <random>

#include <utility>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/generators/EdgeSwitchingMarkovChainGenerator.hpp>
#include <networkit/generators/LFRGenerator.hpp>
#include <networkit/generators/PowerlawDegreeSequence.hpp>
#include <networkit/generators/PubWebGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>

NetworKit::LFRGenerator::LFRGenerator(NetworKit::count n) :
n(n), hasDegreeSequence(false), hasCommunitySizeSequence(false), hasInternalDegreeSequence(false), hasGraph(false), hasPartition(false) { }

void NetworKit::LFRGenerator::setDegreeSequence(std::vector< NetworKit::count > degreeSequence) {
    if (degreeSequence.size() != n) throw std::runtime_error("The degree sequence must have as many entries as there are nodes");
    if (*std::max_element(degreeSequence.begin(), degreeSequence.end()) >= n) throw std::runtime_error("The maximum degree must be smaller than the number of nodes");

    this->degreeSequence = std::move(degreeSequence);
    this->hasDegreeSequence = true;
}

void NetworKit::LFRGenerator::generatePowerlawDegreeSequence(NetworKit::count avgDegree, NetworKit::count maxDegree, double nodeDegreeExp) {
    if (maxDegree >= n) throw std::runtime_error("The maximum degree must be smaller than the number of nodes");

    PowerlawDegreeSequence nodeDegreeSequence(1, maxDegree, nodeDegreeExp);

    nodeDegreeSequence.setMinimumFromAverageDegree(avgDegree);
    nodeDegreeSequence.run();

    this->degreeSequence = nodeDegreeSequence.getDegreeSequence(n);
    this->hasDegreeSequence = true;
}

void NetworKit::LFRGenerator::setCommunitySizeSequence(std::vector< NetworKit::count > communitySizeSequence) {
    this->communitySizeSequence = std::move(communitySizeSequence);
    this->hasCommunitySizeSequence = true;
    this->hasPartition = false;
}

void NetworKit::LFRGenerator::generatePowerlawCommunitySizeSequence(NetworKit::count minCommunitySize, NetworKit::count maxCommunitySize, double communitySizeExp) {
    PowerlawDegreeSequence communityDegreeSequenceGen(minCommunitySize, maxCommunitySize, communitySizeExp);
    communityDegreeSequenceGen.run();

    count sumCommunitySizes = 0;
    communitySizeSequence.clear();

    while (true) {
        count newSize = communityDegreeSequenceGen.getDegree();
        if (sumCommunitySizes + newSize <= n) {
            communitySizeSequence.push_back(newSize);
            sumCommunitySizes += newSize;
        } else { // if the new community doesn't fit anymore, increase the smallest community to fill the gap and exit the loop
            *std::min_element(communitySizeSequence.begin(), communitySizeSequence.end()) += n - sumCommunitySizes;

            break;
        }
    }

    hasCommunitySizeSequence = true;
    this->hasPartition = false;
}

void NetworKit::LFRGenerator::setMu(double mu) {
    if (!hasDegreeSequence) throw std::runtime_error("Error, the degree sequence needs to be set first");

    internalDegreeSequence.resize(n);

    #pragma omp parallel for
    for (omp_index u = 0; u < static_cast<omp_index>(n); ++u) {
        if (degreeSequence[u] == 0) continue;

        double intDeg = (1.0 - mu) * degreeSequence[u];
        if (intDeg < 1) { // assure minimum internal degree of 1
            internalDegreeSequence[u] = 1;
        } else if (Aux::Random::probability() >= std::remainder(intDeg, 1)) {
            internalDegreeSequence[u] = count(intDeg);
        } else {
            internalDegreeSequence[u] = std::ceil(intDeg);
        }
    }

    hasInternalDegreeSequence = true;
}

void NetworKit::LFRGenerator::setMu(const std::vector< double > &mu) {
    if (!hasDegreeSequence) throw std::runtime_error("Error, the degree sequence needs to be set first");
    if (mu.size() != n) throw std::runtime_error("Error, mu must have as many entries as there are nodes");

    internalDegreeSequence.resize(n);

    #pragma omp parallel for
    for (omp_index u = 0; u < static_cast<omp_index>(n); ++u) {
        if (degreeSequence[u] == 0) continue;

        double intDeg = (1.0 - mu[u]) * degreeSequence[u];
        internalDegreeSequence[u] = std::llround(intDeg);
    }

    hasInternalDegreeSequence = true;
}

void NetworKit::LFRGenerator::setMuWithBinomialDistribution(double mu) {
    if (!hasDegreeSequence) throw std::runtime_error("Error, the degree sequence needs to be set first");

    internalDegreeSequence.resize(n);

    auto &gen = Aux::Random::getURNG();
    std::binomial_distribution<count> binDist;

    for (node u = 0; u < n; ++u) {
        if (degreeSequence[u] == 0) continue;

        internalDegreeSequence[u] = binDist(gen, std::binomial_distribution<count>::param_type(degreeSequence[u], 1.0 - mu));
    }

    hasInternalDegreeSequence = true;
}


void NetworKit::LFRGenerator::setPartition(NetworKit::Partition zeta) {
    this->zeta = std::move(zeta);
    this->hasPartition = true;
    this->hasCommunitySizeSequence = false;
}

NetworKit::Graph NetworKit::LFRGenerator::generateIntraClusterGraph(std::vector< NetworKit::count > intraDegreeSequence, const std::vector<NetworKit::node> &localToGlobalNode) {
    // check if sum of degrees is even and fix if necessary
    count intraDegSum = std::accumulate(intraDegreeSequence.begin(), intraDegreeSequence.end(), count{0});

    DEBUG("Possibly correcting the degree sequence");

    for (index j = 0; intraDegSum % 2 != 0 && j < intraDegreeSequence.size(); ++j) {
        index i = Aux::Random::index(intraDegreeSequence.size());
        node u = localToGlobalNode[i];
        if (Aux::Random::real() >= 0.5) {
            if (intraDegreeSequence[i] < intraDegreeSequence.size() - 1 && intraDegreeSequence[i] < degreeSequence[u]) {
                TRACE("Making degree distribution even by increasing intra-degree of node ", u);
                ++intraDegreeSequence[i];
                ++internalDegreeSequence[u];
                ++intraDegSum;
            }
        } else {
            if (intraDegreeSequence[i] > 1) {
                TRACE("Making degree distribution even by decreasing intra-degree of node ", u);
                --intraDegreeSequence[i];
                --internalDegreeSequence[u];
                --intraDegSum;
            }
        }
    }

    DEBUG("Generating intra-cluster graph");
        EdgeSwitchingMarkovChainGenerator intraGen(intraDegreeSequence, true);
    /* even though the sum is even the degree distribution isn't necessarily realizable.
    Disabling the check means that some edges might not be created because of this but at least we will get a graph. */
    return  intraGen.generate();
}


NetworKit::Graph NetworKit::LFRGenerator::generateInterClusterGraph(const std::vector< NetworKit::count > &externalDegreeSequence) {
    EdgeSwitchingMarkovChainGenerator graphGen(externalDegreeSequence, true);

    Graph interG = graphGen.generate();

    // rewire intra-cluster edges in interG as only inter-cluster edges should be generated
    std::vector<std::pair<node, node> > edgesToRemove;

    interG.forEdges([&](node u, node v) {
        if (zeta[u] == zeta[v]) {
            edgesToRemove.emplace_back(u, v);
        }
    });

    if (!edgesToRemove.empty()) {
        std::vector<node> nodeSelection;
        nodeSelection.reserve(interG.numberOfEdges() * 2);

        interG.forNodes([&](node u) {
            for (count i = 0; i < interG.degree(u); ++i) {
                nodeSelection.push_back(u);
            }
        });


        count intraRemovalAttempts = 0;
        count maxIntraRemovelAttempts = interG.numberOfEdges() * 10;
        while (!edgesToRemove.empty()) {
            index i = Aux::Random::index(edgesToRemove.size());
            node s1, t1;
            std::tie(s1, t1) = edgesToRemove[i];

            if (!interG.hasEdge(s1, t1)) {
                edgesToRemove[i] = edgesToRemove.back();
                edgesToRemove.pop_back();
                continue;
            }

            ++intraRemovalAttempts;

            node s2 = Aux::Random::choice(nodeSelection);

            if (s2 == s1 || s2 == t1) continue;

            node t2 = GraphTools::randomNeighbor(interG, s2);

            if (t2 == none) continue;

            if (t1 == t2 || s1 == t2) continue;

            if (interG.hasEdge(s1, t2) || interG.hasEdge(s2, t1)) continue; // FIXME make efficient

            interG.swapEdge(s1, t1, s2, t2);

            edgesToRemove[i] = edgesToRemove.back();
            edgesToRemove.pop_back();

            if (zeta[s1] == zeta[t2]) {
                edgesToRemove.emplace_back(s1, t2);
            }

            if (zeta[s2] == zeta[t1]) {
                edgesToRemove.emplace_back(s2, t1);
            }

            // if t1-t2 is in edgesToRemove, it will be removed in a later iteration of the loop

            if (intraRemovalAttempts > maxIntraRemovelAttempts) {
                break;
            }
        }

        for (auto it = edgesToRemove.begin(); it != edgesToRemove.end(); ++it) {
            if (!interG.hasEdge(it->first, it->second)) {
                *it = edgesToRemove.back();
                edgesToRemove.pop_back();
                if (it == edgesToRemove.end()) break;
            }
        }

        if (!edgesToRemove.empty()) {
            WARN("There are ", edgesToRemove.size(), " intra-cluster edges (which might be counted twice) that should actually be rewired to be inter-cluster edges but couldn't be rewired after ",
                maxIntraRemovelAttempts, " attempts. They will be simply dropped now.");

            for (auto e : edgesToRemove) {
                if (interG.hasEdge(e.first, e.second)) { // these edges might be duplicates/in both directions
                    interG.removeEdge(e.first, e.second);
                }
            }
        }
    }

    return interG;
}

std::vector<std::vector<NetworKit::node>> NetworKit::LFRGenerator::assignNodesToCommunities() {
    std::vector<std::vector<NetworKit::node>> communityNodeList;

    bool assignmentSucceeded = true;
    do {
        assignmentSucceeded = true;
        communityNodeList.clear();
        communityNodeList.resize(communitySizeSequence.size());

        std::vector<count> communitySelection;
        { // generate a random permutation of size(c) copies of each community id c
            communitySelection.reserve(n);

            for (index i = 0; i < communitySizeSequence.size(); ++i) {
                for (node u = 0; u < communitySizeSequence[i]; ++u) {
                    communitySelection.push_back(i);
                }
            }

            std::shuffle(communitySelection.begin(), communitySelection.end(), Aux::Random::getURNG());
        }

        // how many nodes are still missing in the community
        std::vector<count> remainingCommunitySizes(communitySizeSequence);

        // nodes that still need to be assigned to a community
        std::vector<node> nodesToAssign;

        // in the first round, assign each node to a randomly chosen community if possible
        for (node u = 0; u < n; ++u) {
            if (communitySizeSequence[communitySelection[u]] > internalDegreeSequence[u]) {
                communityNodeList[communitySelection[u]].push_back(u);
                --remainingCommunitySizes[communitySelection[u]];
            } else {
                nodesToAssign.push_back(u);
            }
        }

        // now assign a random unassigned node to a random feasible community (chosen at random weighted by community size)
        // if the community is full, remove a node from it
        count attempts = 0;
        count totalAttempts = 0;
        while (!nodesToAssign.empty()) {
            ++totalAttempts;
            index c = Aux::Random::choice(communitySelection);

            index i = Aux::Random::index(nodesToAssign.size());
            node u = nodesToAssign[i];

            { // remove index u from nodesToAssign
                nodesToAssign[i] = nodesToAssign.back();
                nodesToAssign.pop_back();
            }

            // pick another community till the node's internal degree fits into the community
            while (internalDegreeSequence[u] >= communitySizeSequence[c]) {
                // we ensured above that there is always a community in which every node will fit
                c = Aux::Random::choice(communitySelection);
            }

            communityNodeList[c].push_back(u);

            // we are lucky, the community is not full yet
            if (remainingCommunitySizes[c] > 0) {
                --remainingCommunitySizes[c];
                attempts = 0;
            } else { // remove a random node from c
                index r = Aux::Random::index(communityNodeList[c].size());
                nodesToAssign.push_back(communityNodeList[c][r]);
                communityNodeList[c][r] = communityNodeList[c].back();
                communityNodeList[c].pop_back();
                ++attempts;
            }

            if (attempts > 3*n) {
                // merge two smallest communities
                WARN("Needed too long to assign nodes to communities. Either there are too many high-degree nodes or the communities are simply too small.");
                WARN("Changing the community sizes by merging the two smallest communities");
                DEBUG(attempts, " nodes were assigned to full communities, in total, ", totalAttempts, " were made");
                auto minIt = std::min_element(communitySizeSequence.begin(), communitySizeSequence.end());
                count minVal = *minIt;

                // delete minimum element
                *minIt = communitySizeSequence.back();
                communitySizeSequence.pop_back();

                minIt = std::min_element(communitySizeSequence.begin(), communitySizeSequence.end());
                *minIt += minVal;

                assignmentSucceeded = false;

                break;
            }
        }
    } while (!assignmentSucceeded);

    return communityNodeList;
}



void NetworKit::LFRGenerator::run() {
    if (!hasDegreeSequence) throw std::runtime_error("Error, the degree sequence needs to be set first");
    if (!(hasCommunitySizeSequence || hasPartition)) throw std::runtime_error("Error, either the community size sequence or the partition needs to be set first");
    if (!hasInternalDegreeSequence) throw std::runtime_error("Error, mu needs to be set first");

    hasRun = false;

    Aux::SignalHandler handler;

    auto minMaxInternalDegreeIt = std::minmax_element(internalDegreeSequence.begin(), internalDegreeSequence.end());
    count minInternalDegree = *minMaxInternalDegreeIt.first, maxInternalDegree = *minMaxInternalDegreeIt.second;

    handler.assureRunning();

    if (hasCommunitySizeSequence && !hasPartition) { // check if community size sequence is realizable
        auto minMaxCommunitySizeIt = std::minmax_element(communitySizeSequence.begin(), communitySizeSequence.end());
        count minCommunitySize = *minMaxCommunitySizeIt.first, maxCommunitySize = *minMaxCommunitySizeIt.second;

        if (maxInternalDegree >= maxCommunitySize) throw std::runtime_error("Graph not realizable, the maximum internal degree is greater than the largest possible internal degree.");
        if (minInternalDegree >= minCommunitySize) throw std::runtime_error("Graph not realizable, no node can be placed in the smallest community.");

        handler.assureRunning();

        std::vector<count> sortedInternalDegreeSequence = internalDegreeSequence;
        Aux::Parallel::sort(sortedInternalDegreeSequence.begin(), sortedInternalDegreeSequence.end());

        handler.assureRunning();

        std::vector<count> sortedCommunitySizes = communitySizeSequence;
        Aux::Parallel::sort(sortedCommunitySizes.begin(), sortedCommunitySizes.end());

        handler.assureRunning();

        auto communitySizeIt = sortedCommunitySizes.begin();
        count nodesInCommunity = 0;
        for (count deg : sortedInternalDegreeSequence) {
            if (nodesInCommunity == *communitySizeIt) {
                ++communitySizeIt;
                nodesInCommunity = 0;
            }

            if (deg >= *communitySizeIt) {
                throw std::runtime_error("Graph not realizable, community sizes too small or internal degrees too large");
            }

            ++nodesInCommunity;
        }
    } else {
        auto communitySizes = zeta.subsetSizeMap();
        for (index u = 0; u < internalDegreeSequence.size(); ++u) {
            if (internalDegreeSequence[u] >= communitySizes[zeta[u]]) {
                throw std::runtime_error("Graph not realizable, internal degree too large for the assigned community");
            }
        }
    }

    hasGraph = false;
    // initialize the result graph
    G = Graph(n);

    handler.assureRunning();

    // a list of nodes for each community
    std::vector<std::vector<node> > communityNodeList;

    // set communityNodeList and zeta (if not given)
    if (hasCommunitySizeSequence && !hasPartition) {
        communityNodeList = assignNodesToCommunities();

        handler.assureRunning();

        // write communityNodeList into partition instance
        zeta = Partition(n);
        zeta.setUpperBound(communityNodeList.size());

        #pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(communityNodeList.size()); ++i) {
            for (node u : communityNodeList[i]) {
                zeta[u] = i;
            }
        }
    } else {
        communityNodeList.resize(zeta.upperBound());
        G.forNodes([&](node u) {
            communityNodeList[zeta[u]].push_back(u);
        });
    }

    // generate intra-cluster edges
    #pragma omp parallel for // note: parallelization only works because the communities are non-overlapping
    for (omp_index i = 0; i < static_cast<omp_index>(communityNodeList.size()); ++i) {
        const auto &communityNodes = communityNodeList[i];
        if (communityNodes.empty()) continue;

        std::vector<count> intraDeg;
        intraDeg.reserve(communityNodes.size());

        for (node u : communityNodes) {
            intraDeg.push_back(internalDegreeSequence[u]);
        }

        handler.assureRunning();

        Graph intraG = generateIntraClusterGraph(std::move(intraDeg), communityNodes);

        handler.assureRunning();

        #pragma omp critical (generators_lfr_intra_to_global)
        { // FIXME: if we used a graph builder here, we would not need any critical section (only needed for global edge counter)
            intraG.forEdges([&](node i, node j) {
                G.addEdge(communityNodes[i], communityNodes[j]);
            });
        }

        handler.assureRunning();
    }

    // generate inter-cluster edges
    std::vector<count> externalDegree(n);

    G.parallelForNodes([&](node u) {
        externalDegree[u] = degreeSequence[u] - internalDegreeSequence[u];
    });

    handler.assureRunning();

    Graph interG = generateInterClusterGraph(externalDegree);

    handler.assureRunning();

    count edgesRemoved = 0;
    interG.forEdges([&](node u, node v) {
#ifndef NDEBUG // check if edge rewiring actually worked.
        if (!G.hasEdge(u, v)) {
#endif
            G.addEdge(u, v);
#ifndef NDEBUG
        } else {
            ++edgesRemoved;
        }
#endif
    });

    handler.assureRunning();

    if (edgesRemoved > 0) {
        ERROR("Removed ", edgesRemoved, " intra-community edges that already existed but should have been rewired actually!");
    }

    hasGraph = true;
    hasPartition = true;
    hasRun = true;
}

NetworKit::Graph NetworKit::LFRGenerator::generate() {
    run();
    return getMoveGraph();
}

NetworKit::Graph NetworKit::LFRGenerator::getGraph() const {
    if (!hasGraph) throw std::runtime_error("Run must be called first");

    return G;
}

NetworKit::Graph && NetworKit::LFRGenerator::getMoveGraph() {
    if (!hasGraph) throw std::runtime_error("Run must be called first");

    hasGraph = false;

    return std::move(G);
}

NetworKit::Partition NetworKit::LFRGenerator::getPartition() const {
    if (!hasPartition) throw std::runtime_error("Run must be called first");

    return zeta;
}

NetworKit::Partition && NetworKit::LFRGenerator::getMovePartition() {
    if (!hasPartition) throw std::runtime_error("Run must be called first");

    return std::move(zeta);
}

std::string NetworKit::LFRGenerator::toString() const {
    std::stringstream stream;
    stream << "LFRGenerator(" << n << ")";
    return stream.str();
}

bool NetworKit::LFRGenerator::isParallel() const {
    return false;
}
