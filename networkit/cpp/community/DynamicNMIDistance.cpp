// no-networkit-format
/*
 * DynamicNMIDistance.cpp
 *
 *  Created on: Jun 26, 2013
 *      Author: Henning
 */

#include <networkit/community/DynamicNMIDistance.hpp>
#include <networkit/auxiliary/MissingMath.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/Log.hpp>

#include <tlx/unused.hpp>

namespace NetworKit {

bool DynamicNMIDistance::isInBoth(node u, const Partition& oldClustering, const Partition& newClustering) {
    return ((newClustering[u] != none) &&
            (u < oldClustering.numberOfElements()) &&
            (oldClustering[u] != none));
            // number of entries that actually does not exist in Clustering.h
}


/**
 * Formula follows Dhillon, Guan, Kulis: A Unified View of Kernel k-means, ...
 */
double DynamicNMIDistance::getDissimilarity(const Graph& newGraph,
        const Partition& oldClustering, const Partition& newClustering) {

    INFO("compressing clusterings");
//	oldClustering.compact();
//	newClustering.compact();
    INFO("calculating dissimilarity");

    auto log_b = Aux::MissingMath::log_b; // import convenient logarithm function

//	count n = newGraph.numberOfNodes();

    DEBUG("oldClustering=" , oldClustering.getVector());
    DEBUG("newClustering=" , newClustering.getVector());


    std::vector<count> size_old(oldClustering.upperBound());
    std::vector<count> size_new(newClustering.upperBound());

    // precompute sizes for each cluster
    newGraph.forNodes([&](node u){
        if (isInBoth(u, oldClustering, newClustering)) {
            index C = oldClustering[u];
            index D = newClustering[u];
            size_old[C]++;
            size_new[D]++;
        }
    });

    DEBUG("size_old=", size_old);
    DEBUG("size_new=", size_new);


    // confusion matrix
    std::vector<std::vector<count> > confMatrix = this->confusionMatrix(newGraph, oldClustering, newClustering);

    auto numOverlap = [&](Matrix& confMatrix) {
        count num = 0;
        for (Matrix::iterator iter = confMatrix.begin(); iter != confMatrix.end(); ++iter) {
            for (index i = 0; i < iter->size(); ++i) {
                num += (*iter)[i];
            }
        }
        return num;
    };

    count totalOverlap = numOverlap(confMatrix);
    double numDouble = (double) totalOverlap;

    double MI = 0.0; // mutual information
    for (index C = 0; C < oldClustering.upperBound(); C++) {
        for (index D = 0; D < newClustering.upperBound(); D++) {
            count currOverlap = confMatrix[C][D];
            if (currOverlap > 0) {
                double factor1 = (double) currOverlap / (double) numDouble;
                double nominator = (double) (currOverlap * numDouble);
                double aggregate1 = (double) size_old[C]; //  compAggregate1(confMatrix, oldClustering, D);
                double aggregate2 = (double) size_new[D]; // compAggregate2(confMatrix, C, newClustering);
                double denom = aggregate1 * aggregate2;
                DEBUG("frac: " , nominator , " / " , denom , " = " , nominator / denom);
                double factor2 = log_b(nominator / denom, 2);
                MI += factor1 * factor2;

                DEBUG("contribution of " , C , " and " , D , ": " , factor1 , " * " , factor2 , " = " , factor1 * factor2);
            }
        }
    }


    // precompute cluster probabilities
    std::vector<double> P_old(oldClustering.upperBound(), 0.0);
    std::vector<double> P_new(newClustering.upperBound(), 0.0);
    #pragma omp parallel for
    for (omp_index C = static_cast<omp_index>(oldClustering.lowerBound());
         C < static_cast<omp_index>(oldClustering.upperBound()); ++C) {
        P_old[C] = ((double) size_old[C]) / numDouble;
    }
    #pragma omp parallel for
    for (omp_index C = static_cast<omp_index>(newClustering.lowerBound());
         C < static_cast<omp_index>(newClustering.upperBound()); ++C) {
        P_new[C] = ((double) size_new[C]) / numDouble;
    }

    // sanity check
    assert (! std::isnan(MI));
    assert (MI >= 0.0);


    // compute entropy for both clusterings
    double H_old = entropy(oldClustering, totalOverlap, P_old);
    double H_new = entropy(newClustering, totalOverlap, P_new);

    // calculate NMID:
    /* $NMI(\zeta,\eta):=\frac{2\cdot MI(\zeta,\eta)}{H(\zeta)+H\text{(\eta)}}$
     * $NMID(\zeta,\eta):=\begin{cases}
     *	1-NMI(\zeta,\eta)\\
     *	0 & H(\zeta)+H(\eta)=0
     *	\end{cases}$$
     */
    double NMID = 0.0;
    double NMI = 0.0;
    double H_sum = H_old + H_new;
    combineValues(H_sum, MI, NMI, NMID);
    sanityCheck(NMI, NMID);

    return NMID;
}

void DynamicNMIDistance::combineValues(double H_sum, double MI, double& NMI, double& NMID) const {
    if (Aux::NumericTools::equal(H_sum, 0.0)) {
        NMID = 0.0;
    } else {
        NMI = (2.0 * MI) / H_sum;
        NMID = 1.0 - NMI;
    }
}

double DynamicNMIDistance::entropy(const Partition& clustering, count n, std::vector<double> probs) {
    auto log_b = Aux::MissingMath::log_b; // import convenient logarithm function

    // $H(\zeta):=-\sum_{C\in\zeta}P(C)\cdot\log_{2}(P(C))$
    double H = 0.0;
    #pragma omp parallel for reduction(+:H)
    for (omp_index C = static_cast<omp_index>(clustering.lowerBound());
         C < static_cast<omp_index>(clustering.upperBound()); ++C) {
        if (probs[C] != 0) {
            H += probs[C] * log_b(probs[C], 2);
        } // log(0) is not defined
    }
    H = -1.0 * H;

    assert (! std::isnan(H));

    // entropy values range from 0 for the 1-clustering to log_2(n) for the singleton clustering
    assert (Aux::NumericTools::ge(H, 0.0));
    assert (Aux::NumericTools::le(H, log_b(n, 2)));
  (void)n;

    return H;
}

void DynamicNMIDistance::sanityCheck(double& NMI, double& NMID) const {
    DEBUG("sanity check, NMI: " , NMI);
    tlx::unused(NMI);

    if (Aux::NumericTools::equal(NMID, 0.0)) {
        NMID = 0.0;
    }
    if (Aux::NumericTools::equal(NMID, 1.0)) {
        NMID = 1.0;
    }

    // if NMID is close to 0 because of numerical error
    if (! Aux::NumericTools::ge(NMID, 0.0)) {
        ERROR("Set NMID from below 0 to exactly 0: " , NMID);
        NMID = 0.0;
    }
    if (! Aux::NumericTools::le(NMID, 1.0)) {
        ERROR("Set NMID larger than 1 to exactly 1: " , NMID);
        NMID = 1.0;
    }

    assert (Aux::NumericTools::ge(NMID, 0.0));
    assert (Aux::NumericTools::le(NMID, 1.0));
}

std::vector<std::vector<count> > DynamicNMIDistance::confusionMatrix(const Graph &,
        const Partition& first, const Partition& second) {
    index firstUpperId = first.upperBound();
    index secondUpperId = second.upperBound();
    std::vector<std::vector<count> > confMatrix(firstUpperId);

    for (index i = 0; i < first.upperBound(); ++i) {
        confMatrix[i].resize(secondUpperId, 0);
    }

    TRACE("upperId in first, second: " , first.upperBound() , ", " , secondUpperId);

    second.forEntries([&](node u, index secondId) {
        if (this->isInBoth(u, first, second)) {
            TRACE("node " , u , ", id in first: " , first[u] , ", in second: " , second[u]);
            index firstId = first[u];
            assert(firstId < confMatrix.size() && secondId < confMatrix[firstId].size());
            confMatrix[firstId][secondId]++;
        }
    });

    return confMatrix;
}



} /* namespace NetworKit */
