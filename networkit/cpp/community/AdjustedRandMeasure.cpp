// no-networkit-format
#include <networkit/community/AdjustedRandMeasure.hpp>
#include <networkit/community/PartitionIntersection.hpp>


double NetworKit::AdjustedRandMeasure::getDissimilarity(const NetworKit::Graph &G, const NetworKit::Partition &zeta, const NetworKit::Partition &eta) {
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


    count randIndex = 0;
    for (count s : size_intersection) {
        randIndex += s * (s - 1) / 2;
    }

    count sumZeta = 0;
    for (count s : size_zeta) {
        sumZeta += s * (s - 1) / 2;
    }

    count sumEta = 0;
    for (count s : size_eta) {
        sumEta += s * (s - 1) / 2;
    }

    const count denominator = (G.numberOfNodes() * (G.numberOfNodes() - 1)) / 2;

    double maxIndex = 0.5 * static_cast<double>(sumZeta + sumEta);

    double expectedIndex = static_cast<double>(sumZeta * sumEta) / static_cast<double>(denominator);

    if (maxIndex == 0) { // both clusterings are singleton clusterings
        return 0.0;
    } else if (maxIndex == expectedIndex) { // both partitions contain one cluster the whole graph
        return 0.0;
    } else {
        return 1.0 - (randIndex - expectedIndex) / (maxIndex - expectedIndex);
    }
}
