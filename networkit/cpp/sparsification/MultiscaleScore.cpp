// no-networkit-format
/*
 * MultiscaleScore.cpp
 *
 *  Created on: 20.06.2014
 *      Author: Gerd Lindner
 */

#include <networkit/sparsification/MultiscaleScore.hpp>

namespace NetworKit {

MultiscaleScore::MultiscaleScore(const Graph& G, const std::vector<double>& attribute) : EdgeScore<double>(G), attribute(&attribute) {}

void MultiscaleScore::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    //The following vector is used for the _local_ normalization of edgeweights.
    //We use a global vector for performance reasons.
    std::vector<edgeweight> normalizedWeights(G->upperNodeIdBound());

    std::vector<double> multiscaleAttribute(G->upperEdgeIdBound(), 0.0);

    G->forNodes([&](node u) {
        count k = G->degree(u);

        //Normalize edgeweights of N(u)
        edgeweight sum = 0.0;
        G->forNeighborsOf(u, [&](node, node, edgeid eid) {
            sum += (*attribute)[eid];
        });
        G->forNeighborsOf(u, [&](node, node v, edgeid eid) {
            normalizedWeights[v] = (*attribute)[eid] / sum;
        });

        //Filter edges by probability
        G->forNeighborsOf(u, [&](node, node v, edgeid eid) {
            //In case d(u) == 1 and d(v) > 1: ignore u
            //if (k > 1 || G.degree(v) == 1) {
                edgeweight p = normalizedWeights[v];
                double probability = getProbability(k, p);

                multiscaleAttribute[eid] = std::max(multiscaleAttribute[eid], probability);
            //}
        });
    });

    scoreData = std::move(multiscaleAttribute);
    hasRun = true;
}

/**
 * Returns the probability that a node of the given
 * degree has an edge of the given weight.
 *
 * The null hypothesis is the following: the normalized weights of the
 * edges connected to a node of degree k are uniformly distributed.
 */
double MultiscaleScore::getProbability(count degree, edgeweight normalizedWeight) {
    return 1.0 - pow(1.0 - normalizedWeight, static_cast<double>(degree) - 1.0);
}

double MultiscaleScore::score(node, node) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

double MultiscaleScore::score(edgeid) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
