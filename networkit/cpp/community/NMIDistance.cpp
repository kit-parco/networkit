// no-networkit-format
/*
 * NormalizedMutualInformation.cpp
 *
 *  Created on: 30.04.2013
 *      Author: cls
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/MissingMath.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/community/DynamicNMIDistance.hpp>
#include <networkit/community/NMIDistance.hpp>
#include <networkit/community/PartitionIntersection.hpp>

namespace NetworKit {

double NMIDistance::getDissimilarity(const Graph& G, const Partition& zeta, const Partition& eta) {

    double n = G.numberOfNodes();

    DEBUG("zeta=" , zeta.getVector());
    DEBUG("eta=" , eta.getVector());

    std::vector<count> size_zeta(zeta.upperBound(), 0), size_eta(eta.upperBound(), 0);

    // precompute sizes for each cluster

    G.forNodes([&](node u){
        index C = zeta[u];
        index D = eta[u];
        assert (C != none);
        assert (D != none);
        size_zeta[C] += 1;
        size_eta[D] += 1;
    });

    DEBUG("size_zeta=" , size_zeta);
    DEBUG("size_eta=" , size_eta);

    // precompute cluster probabilities
    std::vector<double> P_zeta(zeta.upperBound(), 0.0);
    std::vector<double> P_eta(eta.upperBound(), 0.0);

    for (index C = 0; C < size_zeta.size(); ++C) {
        P_zeta[C] = static_cast<double>(size_zeta[C]) / n;
    }

    for (index D = 0; D < size_eta.size(); ++D) {
        P_eta[D] = static_cast<double>(size_eta[D]) / n;
    }

    Partition overlap = PartitionIntersection().calculate(zeta, eta);
    DEBUG("overlap=", overlap.getVector());

    std::vector<index> overlap_zeta(overlap.upperBound(), none), overlap_eta(overlap.upperBound(), none); // map from the overlap to zeta and eta

    std::vector<count> overlapSizes(overlap.upperBound(), 0); // overlapSizes[O] returns the size of the overlap cluster

    G.forNodes([&](node u){
        index O = overlap[u];
        if (overlapSizes[O] == 0) {
            overlap_zeta[O] = zeta[u];
            overlap_eta[O] = eta[u];
        }
        overlapSizes[O] += 1;
    });

    DEBUG("overlapSizes=", overlapSizes);


    auto log_b = Aux::MissingMath::log_b; // import convenient logarithm function


    // calculate mutual information
    // $MI(\zeta,\eta):=\sum_{C\in\zeta}\sum_{D\in\eta}\frac{|C\cap D|}{n}\cdot\log_{2}\left(\frac{|C\cap D|\cdot n}{|C|\cdot|D|}\right)$
    double MI = 0.0; // mutual information
    for (index O = 0; O < overlapSizes.size(); ++O) {
        if (overlapSizes[O] > 0) {
            index C = overlap_zeta[O];
            index D = overlap_eta[O];
            count sizeC = size_zeta[C];
            count sizeD = size_eta[D];
            count sizeO = overlapSizes[O];
            double factor1 =  static_cast<double>(sizeO) / n;
            assert ((sizeC * sizeD) != 0);
            TRACE("overlap of " , C , " and " , D , " has size: " , sizeO);
            TRACE("union of " , C , " and " , D , " has size: " , (sizeD + sizeC - sizeO));
            double frac2 = (static_cast<double>(sizeO) * n) / static_cast<double>(sizeC * sizeD);
            assert (frac2 != 0);
            double factor2 = log_b(frac2, 2);
            TRACE("frac2 = " , frac2 , ", factor1 = " , factor1 , ", factor2 = " , factor2);
            MI += factor1 * factor2;
        }
    }

    // sanity check
    assert (! std::isnan(MI));
    assert (MI >= 0.0);

    // compute entropy for both clusterings
    DynamicNMIDistance dynNMID;
    double H_zeta = dynNMID.entropy(zeta, n, P_zeta);
    double H_eta = dynNMID.entropy(eta, n, P_eta);

    assert (! std::isnan(H_zeta));
    assert (! std::isnan(H_eta));

    if (H_zeta == 0 || H_eta == 0) {
        WARN("You are comparing two partitions where one has zero entropy (i.e. consists of one cluster), ",
        "thus the mutual information will always be zero which might not be what you intended. ",
        "Note that two partitions will still be reported to be equal when both have zero mutual information.");
    }

    double H_sum = H_zeta + H_eta;
    double NMI = 0.0;
    double NMID = 0.0;
    dynNMID.combineValues(H_sum, MI, NMI, NMID);
    dynNMID.sanityCheck(NMI, NMID);

    return NMID;
}

} /* namespace NetworKit */
