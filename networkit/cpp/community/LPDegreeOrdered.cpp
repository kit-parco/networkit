// no-networkit-format
/*
 * LPDegreeOrdered.cpp
 *
 *  Created on: 24.09.2013
 *      Author: cls
 */

#include <unordered_map>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/community/LPDegreeOrdered.hpp>

namespace NetworKit {

LPDegreeOrdered::LPDegreeOrdered(const Graph& G) : CommunityDetectionAlgorithm(G) {}

void LPDegreeOrdered::run() {
    count n = G->numberOfNodes();
    count theta = n / 1e5;
    DEBUG("theta: " , theta);

    index z = G->upperNodeIdBound();
    Partition labels(z);
    // initialize all labels to singletons
    labels.allToSingletons();

    // initialize all nodes as active
    std::vector<int> active(z + 1, 1); // not a boolean vector because there might be problems with parallel access

    count nUpdated; // number of nodes which have been updated in last iteration
    nUpdated = n; // all nodes have new labels -> first loop iteration runs
    nIterations = 0; // number of iterations

    auto propagateLabels = [&](node v){
        if ((active[v]) && (G->degree(v) > 0)) {
            std::unordered_map<label, count> labelCounts; // neighborLabelCounts maps label -> frequency in the neighbors
            // count the labels in the neighborhood of v
            G->forNeighborsOf(v, [&](node w) {
                label lw = labels.subsetOf(w);//labels[w];
                labelCounts[lw] += 1; // add weight of edge {v, w}
            });

            // get dominant label
            label dominant = std::max_element(labelCounts.begin(),
                            labelCounts.end(),
                            [](const std::pair<label, count>& p1, const std::pair<label, count>& p2) {
                                return p1.second < p2.second;})->first;
            if (labels[v] != dominant) { // UPDATE
                labels.moveToSubset(dominant,v);//labels[v] = dominant;
                nUpdated += 1; // TODO: atomic update?
                G->forNeighborsOf(v, [&](node u) {
                    active[u] = 1;
                });
            } else {
                active[v] = 0;
            }

        } // else node is isolated or inactive
    };


    // sort nodes by degree
    std::vector<node> nodes;
    G->forNodes([&](node v) {
        nodes.push_back(v);
    });

    Aux::Parallel::sort(nodes.begin(), nodes.end(), [&](node u, node v) {
        return G->degree(u) < G->degree(v); // lower degree before higher degree
    });

    // propagate labels
    while (nUpdated > theta) { // as long as a label has changed...
        nUpdated = 0; // reset update counter

        for (node v : nodes) {
            propagateLabels(v);
        }

        INFO("updated labels: " , nUpdated);
        nIterations += 1;
    }
    result = std::move(labels);
    hasRun = true;
}


count LPDegreeOrdered::numberOfIterations() {
    return this->nIterations;
}

std::string LPDegreeOrdered::toString() const {
    return "LPDegreeOrdered()";
}

} /* namespace NetworKit */
