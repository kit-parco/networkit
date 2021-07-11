// no-networkit-format
/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#include <networkit/components/ParallelConnectedComponents.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/auxiliary/Log.hpp>

namespace NetworKit {

ParallelConnectedComponents::ParallelConnectedComponents(const Graph& G, bool coarsening) : Algorithm(), G(&G), coarsening(coarsening) {

}


void ParallelConnectedComponents::run() {
    if (G->isDirected()) {
        throw std::runtime_error("algorithm does not accept directed graphs");
    }

    // calculate connected components by label propagation
    count z = G->upperNodeIdBound();

    DEBUG("initializing labels");
    component = Partition(z);
    component.allToSingletons();
    // remove nodes that do not exist from the partition so it doesn't report wrong numbers
    component.parallelForEntries([&](node u, index) {
        if (!G->hasNode(u)) {
            component[u] = none;
        }
    });

    DEBUG("initializing active nodes");
    const char INACTIVE = 0;
    const char ACTIVE = 1;
    std::vector<char> activeNodes(z); // record if node must be processed
    std::vector<char> nextActiveNodes(z, ACTIVE); // for next iteration

    DEBUG("main loop");
//	count numActive = 0; // for debugging purposes only
    count numIterations = 0;
    bool change = true;
    // only 8 iterations when coarsening is on, otherwise till no more changes happened
    while (change && (!coarsening || numIterations < 8)) {
//		TRACE("label propagation iteration");
        activeNodes.swap(nextActiveNodes);
        nextActiveNodes.assign(z, INACTIVE);
        change = false;
//		numActive = 0;
        G->balancedParallelForNodes([&](node u) {
            if (activeNodes[u] == ACTIVE) {
//				++numActive;
                // get smallest
                index smallest = component[u];
                G->forNeighborsOf(u, [&](node v) {
                    smallest = std::min(smallest, component[v]);
                });

                if (component[u] != smallest) {
                    component.moveToSubset(smallest, u);
                    change = true;
                    G->forNeighborsOf(u, [&](node v) {
                        // only nodes that do not have the smallest component label
                        // will see a new component label because of the change of u,
                        // only they need to be activated
                        if (component[v] != smallest) {
                            nextActiveNodes[v] = ACTIVE;
                        }
                    });
                }
            }
        });
//		TRACE("num active: ", numActive);
        ++numIterations;
    }
    if (coarsening && numIterations == 8) { // TODO: externalize constant
        // coarsen and make recursive call
        ParallelPartitionCoarsening con(*G, component);
        con.run();
        auto Gcon = con.getCoarseGraph();
        ParallelConnectedComponents cc(Gcon);
        cc.run();

        // apply to current graph
        auto nodeMapping = con.getFineToCoarseNodeMapping();
        G->parallelForNodes([&](node u) {
            component[u] = cc.componentOfNode(nodeMapping[u]);
        });
    }

    hasRun = true;
}

void ParallelConnectedComponents::runSequential() {
    if (G->isDirected()) {
        throw std::runtime_error("algorithm does not accept directed graphs");
    }
    // calculate connected components by label propagation
    count z = G->upperNodeIdBound();
    DEBUG("initializing labels");
    component = Partition(z);
    component.allToSingletons();
    // remove nodes that do not exist from the partition so it doesn't report wrong numbers
    component.forEntries([&](node u, index) {
        if (!G->hasNode(u)) {
            component[u] = none;
        }
    });

    DEBUG("initializing active nodes");
    std::vector<bool> activeNodes(z, true); // record if node must be processed

    DEBUG("main loop");
    // count numActive = 0; // for debugging purposes only
    count numIterations = 0;
    bool change = true;
    // only 8 iterations when coarsening is on, otherwise till no more changes happened
    while (change && (!coarsening || numIterations < 8)) {
        // TRACE("label propagation iteration");
        change = false;
        // numActive = 0;
        G->forNodes([&](node u) {
            if (activeNodes[u]) {
                // ++numActive;
                // get smallest
                index smallest = component[u];
                G->forNeighborsOf(u, [&](node v) {
                    smallest = std::min(smallest, component[v]);
                });

                if (component[u] != smallest) {
                    component.moveToSubset(smallest, u);
                    change = true;
                    G->forNeighborsOf(u, [&](node v) {
                        // only nodes that do not have the smallest component label
                        // will see a new component label because of the change of u,
                        // only they need to be activated
                        if (component[v] != smallest) {
                            activeNodes[v] = true;
                        }
                    });
                } else {
                    activeNodes[u] = false; // current node becomes inactive
                }
            }
        });
        // TRACE("num active: ", numActive);
        ++numIterations;
    }
    if (coarsening && numIterations == 8) { // TODO: externalize constant
        // coarsen and make recursive call
        ParallelPartitionCoarsening con(*G, component);
        con.run();
        auto Gcon = con.getCoarseGraph();
        ParallelConnectedComponents cc(Gcon);
        cc.run();

        // apply to current graph
        auto nodeMapping = con.getFineToCoarseNodeMapping();
        G->forNodes([&](node u) {
            component[u] = cc.componentOfNode(nodeMapping[u]);
        });
    }

    hasRun = true;
}


Partition ParallelConnectedComponents::getPartition() const {
    assureFinished();
    return this->component;
}

count ParallelConnectedComponents::numberOfComponents() const {
    assureFinished();
    return this->component.numberOfSubsets();
}

count ParallelConnectedComponents::componentOfNode(node u) const {
    assureFinished();
    assert (component[u] != none);
    return component[u];
}

std::vector<std::vector<node>> ParallelConnectedComponents::getComponents() const {

    std::vector<std::vector<node>> result(this->numberOfComponents());
    G->forNodes([&](node u) { result[component[u]].push_back(u); });
    return result;
}

}
