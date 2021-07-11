// no-networkit-format
/*
*  EffectiveDiameterApproximation.cpp
*
*  Created on: 29.03.16
*      Author: Maximilian Vogel
*/

#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/EffectiveDiameterApproximation.hpp>

namespace NetworKit {
EffectiveDiameterApproximation::EffectiveDiameterApproximation(const Graph& G, const double ratio, const count k, const count r) : Algorithm(), G(&G), ratio(ratio), k(k), r(r)  {
    if (G.isDirected()) throw std::runtime_error("current implementation can only deal with undirected graphs");
    ConnectedComponents cc(G);
    cc.run();
    if (cc.getPartition().numberOfSubsets() > 1) throw std::runtime_error("current implementation only runs on graphs with 1 connected component");
}

void EffectiveDiameterApproximation::run() {
    count z = G->upperNodeIdBound();
    // the length of the bitmask where the number of connected nodes is saved
    count lengthOfBitmask = (count) ceil(log2(G->numberOfNodes()));
    // saves all k bitmasks for every node of the current iteration
    std::vector<std::vector<unsigned int>> mCurr(z);
    // saves all k bitmasks for every node of the previous iteration
    std::vector<std::vector<unsigned int>> mPrev(z);
    // the maximum possible bitmask based on the random initialization of all k bitmasks
    std::vector<count> highestCount;
    // the amount of nodes that need to be connected to all others nodes
    count threshold = (count) (ceil(ratio * G->numberOfNodes()));
    // the current distance of the neighborhoods
    count h = 1;
    // sums over the number of edges needed to reach 90% of all other nodes
    effectiveDiameter = 0;
    // the estimated number of connected nodes
    double estimatedConnectedNodes;
    // used for setting a random bit in the bitmasks
    double random;
    // the position of the bit that has been set in a bitmask
    count position;
    // nodes that are not connected to enough nodes yet
    std::vector<node> activeNodes;

    // initialize all vectors
    highestCount.assign(k, 0);
    G->forNodes([&](node v) {
        std::vector<unsigned int> bitmasks;
        bitmasks.assign(k, 0);
        mCurr[v] = bitmasks;
        mPrev[v] = bitmasks;
        activeNodes.push_back(v);
        // set one bit in each bitmask with probability P(bit i=1) = 0.5^(i+1), i=0,..
        for (count j = 0; j < k; j++) {
            random = Aux::Random::real(0,1);
            position = ceil(log(random)/log(0.5) - 1);
            // set the bit in the bitmask
            if (position < lengthOfBitmask+r) {
                mPrev[v][j] |= 1 << position;
            }
            // add the current bit to the maximum-bitmask
            highestCount[j] = highestCount[j] | mPrev[v][j];
        }
    });

    // as long as we need to connect more nodes
    while (!activeNodes.empty()) {
        for (count x = 0; x < activeNodes.size(); ++x) {
            node v = activeNodes[x];
            #pragma omp parallel for
            // for each parallel approximation
            for (omp_index j = 0; j < k; j++) {
                // the node is still connected to all previous neighbors
                mCurr[v][j] = mPrev[v][j];
                // and to all previous neighbors of all its neighbors
                G->forNeighborsOf(v, [&](node u) {
                    mCurr[v][j] = mCurr[v][j] | mPrev[u][j];
                });
            }

            // the least bit number in the bitmask of the current node/distance that has not been set
            double b = 0;

            for (count j = 0; j < k; j++) {
                for (count i = 0; i < sizeof(i)*8; i++) {
                    if (((mCurr[v][j] >> i) & 1) == 0) {
                        b += i;
                        break;
                    }
                }
            }
            // calculate the average least bit number that has not been set over all parallel approximations
            b = b / k;

            // calculate the estimated number of neighbors where 0.77351 is a correction factor and the result of a complex sum
            estimatedConnectedNodes = (pow(2,b) / 0.77351);

            // check whether all k bitmask for this node have reached their highest possible value
            bool nodeFinished = true;
            for (count j = 0; j < k; j++) {
                if (mCurr[v][j] != highestCount[j]) {
                    nodeFinished = false;
                    break;
                }
            }
            // if the node wont change or is connected to enough nodes it must no longer be considered
            if (estimatedConnectedNodes >= threshold || nodeFinished) {
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

double EffectiveDiameterApproximation::getEffectiveDiameter() const {
    if(!hasRun) {
        throw std::runtime_error("Call run()-function first.");
    }
    return effectiveDiameter;
}

} /* namespace NetworKit */
