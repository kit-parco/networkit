// no-networkit-format
/*
*  EffectiveDiameter.cpp
*
*  Created on: 16.06.2014
*      Author: Marc Nemes
*/

#include <cmath>
#include <omp.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/EffectiveDiameter.hpp>

namespace NetworKit {

EffectiveDiameter::EffectiveDiameter(const Graph& G, const double ratio) : Algorithm(), G(&G), ratio(ratio) {
    if (G.isDirected()) throw std::runtime_error("current implementation can only deal with undirected graphs");
    ConnectedComponents cc(G);
    cc.run();
    if (cc.numberOfComponents() > 1) throw std::runtime_error("current implementation only runs on graphs with 1 connected component");
}

void EffectiveDiameter::run() {
    count z = G->upperNodeIdBound();
    // saves the reachable nodes of the current iteration
    std::vector<std::vector<bool>> mCurr(z);
    // saves the reachable nodes of the previous iteration
    std::vector<std::vector<bool>> mPrev(z);
    // sums over the number of edges needed to reach 90% of all other nodes
    effectiveDiameter = 0;
    // the current distance of the neighborhoods
    count h = 1;
    // number of nodes that need to be connected with all other nodes
    auto threshold = static_cast<count>(std::ceil(ratio * G->numberOfNodes()));
    // nodes that are not connected to enough nodes yet
    std::vector<node> activeNodes;

    // initialize all nodes
    G->forNodes([&](node v){
        std::vector<bool> connectedNodes;
        // initialize n entries with value 0
        connectedNodes.assign(z, 0);
        // the node is always connected to itself
        connectedNodes[v] = 1;
        mCurr[v] = connectedNodes;
        mPrev[v] = connectedNodes;
        activeNodes.push_back(v);
    });

    // as long as we need to connect more nodes
    while (!activeNodes.empty()) {
        for (count x = 0; x < activeNodes.size(); x++) {
            node v = activeNodes[x];
                mCurr[v] = mPrev[v];
                G->forNeighborsOf(v, [&](node u) {
                    for (count i = 0; i < G->numberOfNodes(); i++) {
                        // add the current neighbor of u to the neighborhood of v
                        mCurr[v][i] = mCurr[v][i] || mPrev[u][i];
                    }
                });

                // compute the number of connected nodes
                count numConnectedNodes = 0;
                for (count i = 0; i < G->numberOfNodes(); i++) {
                    if (mCurr[v][i] == 1) {
                        numConnectedNodes++;
                    }
                }

                // when the number of connected nodes surpasses the threshold the node must no longer be considered
                if (numConnectedNodes >= threshold) {
                    effectiveDiameter += h;
                    // remove the current node from future iterations
                    std::swap(activeNodes[x], activeNodes.back());
                    activeNodes.pop_back();
                    --x; //don't skip former activeNodes.back() that has been switched to activeNodes[x]
                }
            }
            mPrev = mCurr;
            h++;
    }
    effectiveDiameter /= G->numberOfNodes();
    hasRun = true;
}

double EffectiveDiameter::getEffectiveDiameter() const {
    if(!hasRun) {
        throw std::runtime_error("Call run()-function first.");
    }
    return effectiveDiameter;
}

} // namespace NetworKit
