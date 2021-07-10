/*
 * MLPLM.cpp
 *
 *  Created on: 20.11.2013
 *      Author: cls
 */

#include <omp.h>
#include <sstream>
#include <utility>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/coarsening/ClusteringProjector.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/community/PLM.hpp>

namespace NetworKit {

PLM::PLM(const Graph &G, bool refine, double gamma, std::string par, count maxIter, bool turbo,
         bool recurse, Partition zeta)
    : CommunityDetectionAlgorithm(G), zeta(zeta), parallelism(std::move(par)), refine(refine), gamma(gamma),
      maxIter(maxIter), turbo(turbo), recurse(recurse) {}

PLM::PLM(const Graph& G, const PLM& other) : CommunityDetectionAlgorithm(G), zeta(other.zeta), parallelism(other.parallelism), refine(other.refine), gamma(other.gamma), maxIter(other.maxIter), turbo(other.turbo), recurse(other.recurse) {}

void PLM::run() {
    Aux::SignalHandler handler;

    count z = G->upperNodeIdBound();

    // init communities to singletons if none given as input
    if (zeta.numberOfSubsets() == 0) {
        zeta = Partition(z);
        zeta.allToSingletons();
    }
    
    index o = zeta.upperBound();


    // init graph-dependent temporaries
    std::vector<double> volNode(z, 0.0);
    // $\omega(E)$
    edgeweight total = G->totalEdgeWeight();
    DEBUG("total edge weight: " , total);
    edgeweight divisor = (2 * total * total); // needed in modularity calculation

    G->parallelForNodes([&](node u) { // calculate and store volume of each node
        volNode[u] += G->weightedDegree(u);
        volNode[u] += G->weight(u, u); // consider self-loop twice
    });

    // init community-dependent temporaries
    std::vector<double> volCommunity(o, 0.0);
    zeta.parallelForEntries([&](node u, index C) { // set volume for all communities
        if (C != none)
            #pragma omp atomic
            volCommunity[C] += volNode[u];
    });

    // first move phase
    bool moved = false; // indicates whether any node has been moved in the last pass
    bool change = false; // indicates whether the communities have changed at all

    // stores the affinity for each neighboring community (index), one vector per thread
    std::vector<std::vector<edgeweight> > turboAffinity;
    // stores the list of neighboring communities, one vector per thread
    std::vector<std::vector<index> > neigh_comm;


    if (turbo) {
        if (this->parallelism != "none" && this->parallelism != "none randomized") { // initialize arrays for all threads only when actually needed
            turboAffinity.resize(omp_get_max_threads());
            neigh_comm.resize(omp_get_max_threads());
            for (auto &it : turboAffinity) {
                // resize to maximum community id
                it.resize(zeta.upperBound());
            }
        } else { // initialize array only for first thread
            turboAffinity.emplace_back(zeta.upperBound());
            neigh_comm.emplace_back();
        }
    }

    // try to improve modularity by moving a node to neighboring clusters
    auto tryMove = [&](node u) {
        // trying to move node u
        index tid = omp_get_thread_num();

        // collect edge weight to neighbor clusters
        std::map<index, edgeweight> affinity;

        if (turbo) {
            neigh_comm[tid].clear();
            G->forNeighborsOf(u, [&](node v) {
                turboAffinity[tid][zeta[v]] = -1; // set all to -1 so we can see when we get to it the first time
            });
            turboAffinity[tid][zeta[u]] = 0;
            G->forNeighborsOf(u, [&](node v, edgeweight weight) {
                if (u != v) {
                    index C = zeta[v];
                    if (turboAffinity[tid][C] == -1) {
                        // found the neighbor for the first time, initialize to 0 and add to list of neighboring communities
                        turboAffinity[tid][C] = 0;
                        neigh_comm[tid].push_back(C);
                    }
                    turboAffinity[tid][C] += weight;
                }
            });
        } else {
            G->forNeighborsOf(u, [&](node v, edgeweight weight) {
                if (u != v) {
                    index C = zeta[v];
                    affinity[C] += weight;
                }
            });
        }

        // sub-functions

        // $\vol(C \ {x})$ - volume of cluster C excluding node x
        auto volCommunityMinusNode = [&](index C, node x) {
            double volC = 0.0;
            double volN = 0.0;
            volC = volCommunity[C];
            if (zeta[x] == C) {
                volN = volNode[x];
                return volC - volN;
            } else {
                return volC;
            }
        };

        auto modGain = [&](node u, index C, index D, edgeweight affinityC, edgeweight affinityD) {
            double volN = 0.0;
            volN = volNode[u];
            double delta = (affinityD - affinityC) / total + this->gamma * ((volCommunityMinusNode(C, u) - volCommunityMinusNode(D, u)) * volN) / divisor;
            return delta;
        };

        index best = none;
        index C = none;
        double deltaBest = -1;

        C = zeta[u];

        if (turbo) {
            edgeweight affinityC = turboAffinity[tid][C];

            for (index D : neigh_comm[tid]) {
                if (D != C) { // consider only nodes in other clusters (and implicitly only nodes other than u)
                    double delta = modGain(u, C, D, affinityC, turboAffinity[tid][D]);

                    if (delta > deltaBest) {
                        deltaBest = delta;
                        best = D;
                    }
                }
            }
        } else {
            edgeweight affinityC = affinity[C];

            for (auto it : affinity) {
                index D = it.first;
                if (D != C) { // consider only nodes in other clusters (and implicitly only nodes other than u)
                    double delta = modGain(u, C, D, affinityC, it.second);
                    if (delta > deltaBest) {
                        deltaBest = delta;
                        best = D;
                    }
                }
            }
        }

        if (deltaBest > 0) { // if modularity improvement possible
            assert (best != C && best != none);// do not "move" to original cluster

            zeta[u] = best; // move to best cluster
            // node u moved

            // mod update
            double volN = 0.0;
            volN = volNode[u];
            // update the volume of the two clusters
            #pragma omp atomic
            volCommunity[C] -= volN;
            #pragma omp atomic
            volCommunity[best] += volN;

            moved = true; // change to clustering has been made
        }
    };

    // performs node moves
    auto movePhase = [&](){
        count iter = 0;
        do {
            moved = false;
            // apply node movement according to parallelization strategy
            if (this->parallelism == "none") {
                G->forNodes(tryMove);
            } else if (this->parallelism == "simple") {
                G->parallelForNodes(tryMove);
            } else if (this->parallelism == "balanced") {
                G->balancedParallelForNodes(tryMove);
            } else if (this->parallelism == "none randomized") {
                G->forNodesInRandomOrder(tryMove);
            } else {
                ERROR("unknown parallelization strategy: " , this->parallelism);
                throw std::runtime_error("unknown parallelization strategy");
            }
            if (moved) change = true;

            if (iter == maxIter) {
                WARN("move phase aborted after ", maxIter, " iterations");
            }
            iter += 1;
        } while (moved && (iter <= maxIter) && handler.isRunning());
        DEBUG("iterations in move phase: ", iter);
    };
    handler.assureRunning();
    // first move phase
    Aux::Timer timer;
    timer.start();

    movePhase();

    timer.stop();
    timing["move"].push_back(timer.elapsedMilliseconds());
    handler.assureRunning();
    if (recurse && change) {
        DEBUG("nodes moved, so begin coarsening and recursive call");

        timer.start();

        std::pair<Graph, std::vector<node>> coarsened = coarsen(*G, zeta); // coarsen graph according to communities

        timer.stop();
        timing["coarsen"].push_back(timer.elapsedMilliseconds());

        PLM onCoarsened(coarsened.first, this->refine, this->gamma, this->parallelism, this->maxIter, this->turbo);
        onCoarsened.run();
        Partition zetaCoarse = onCoarsened.getPartition();

        // get timings
        auto tim = onCoarsened.getTiming();
        for (count t : tim["move"]) {
            timing["move"].push_back(t);
        }
        for (count t : tim["coarsen"]) {
            timing["coarsen"].push_back(t);
        }
        for (count t : tim["refine"]) {
            timing["refine"].push_back(t);
        }

        DEBUG("coarse graph has ", coarsened.first.numberOfNodes(), " nodes and ", coarsened.first.numberOfEdges(), " edges");
        zeta = prolong(coarsened.first, zetaCoarse, *G, coarsened.second); // unpack communities in coarse graph onto fine graph
        // refinement phase
        if (refine) {
            DEBUG("refinement phase");
            // reinit community-dependent temporaries
            o = zeta.upperBound();
            volCommunity.clear();
            volCommunity.resize(o, 0.0);
            zeta.parallelForEntries([&](node u, index C) { // set volume for all communities
                if (C != none) {
                    edgeweight volN = volNode[u];
                    #pragma omp atomic
                    volCommunity[C] += volN;
                }
            });
            // second move phase
            timer.start();

            movePhase();

            timer.stop();
            timing["refine"].push_back(timer.elapsedMilliseconds());

        }
    }
    result = std::move(zeta);
    hasRun = true;
}

std::string PLM::toString() const {
    std::stringstream stream;
    stream << "PLM(";
    stream << parallelism;
    if (refine) {
        stream << "," << "refine";
    }
    stream << "," << "pc";
    if (turbo) {
        stream << "," << "turbo";
    }
    if (!recurse) {
        stream << "," << "non-recursive";
    }
    stream << ")";

    return stream.str();
}

std::pair<Graph, std::vector<node> > PLM::coarsen(const Graph& G, const Partition& zeta) {
    ParallelPartitionCoarsening parCoarsening(G, zeta);
    parCoarsening.run();
    return {parCoarsening.getCoarseGraph(),parCoarsening.getFineToCoarseNodeMapping()};
}

Partition PLM::prolong(const Graph &, const Partition& zetaCoarse, const Graph& Gfine, std::vector<node> nodeToMetaNode) {
    Partition zetaFine(Gfine.upperNodeIdBound());
    zetaFine.setUpperBound(zetaCoarse.upperBound());

    Gfine.forNodes([&](node v) {
        node mv = nodeToMetaNode[v];
        index cv = zetaCoarse[mv];
        zetaFine[v] = cv;
    });

    return zetaFine;
}

std::map<std::string, std::vector<count> > PLM::getTiming() {
    return timing;
}

} /* namespace NetworKit */
