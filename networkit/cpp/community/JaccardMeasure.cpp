// no-networkit-format
/*
 * JaccardMeasure.cpp
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt
 */

#include <networkit/community/JaccardMeasure.hpp>
#include <networkit/community/PartitionIntersection.hpp>

namespace NetworKit {

double JaccardMeasure::getDissimilarity(const Graph& G, const Partition& zeta,
        const Partition& eta) {

    Partition intersection = PartitionIntersection().calculate(zeta, eta);

    std::vector<count> size_zeta(zeta.upperBound(), 0);
    std::vector<count> size_eta(eta.upperBound(), 0);
    std::vector<count> size_intersection(intersection.upperBound(), 0);

    // precompute sizes for each cluster
    G.forNodes([&](node u){
        index C = zeta[u];
        index D = eta[u];
        index I = intersection[u];
        assert (C != none);
        assert (D != none);
        assert (I != none);
        size_zeta[C] += 1;
        size_eta[D] += 1;
        size_intersection[I] += 1;
    });


    count sumIntersection = 0;
    for (count s : size_intersection) {
        sumIntersection += s * (s - 1) / 2;
    }

    count sumZeta = 0;
    for (count s : size_zeta) {
        sumZeta += s * (s - 1) / 2;
    }

    count sumEta = 0;
    for (count s : size_eta) {
        sumEta += s * (s - 1) / 2;
    }

    double n = G.numberOfNodes();

    count s11 = sumIntersection; // number of node pairs for which clusterings aggree
    double s00 = n * (n-1) / 2 + sumIntersection - static_cast<double>(sumZeta + sumEta); // number of node pairs for which clusterings disagree

    double jaccard;
    double divisor = n * (n - 1) - 2.0 * s00;
    if (divisor > 0) {
        jaccard = 1.0 - ((2.0 * s11) / divisor);
    } else {
        jaccard = 0;
    }

    // assert range [0, 1]
    assert (jaccard <= 1.0);
    assert (jaccard >= 0.0);
    return jaccard;

}

} /* namespace NetworKit */
