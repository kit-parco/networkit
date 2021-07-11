// no-networkit-format
/*
 * NodeStructuralRandMeasure.cpp
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt
 */

#include <networkit/community/NodeStructuralRandMeasure.hpp>
#include <networkit/community/PartitionIntersection.hpp>

namespace NetworKit {

double NodeStructuralRandMeasure::getDissimilarity(const Graph& G, const Partition& zeta, const Partition& eta) {
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

    double A = n * (n-1.0) / 2.0 + 2.0 * static_cast<double>(sumIntersection) - static_cast<double>(sumZeta + sumEta);

    double rand = 1 - ((2 * A) * 1.0 / (n * (n-1)));

    // assert range [0, 1]
    assert (rand <= 1.0);
    assert (rand >= 0.0);
    return rand;
}

} /* namespace NetworKit */
