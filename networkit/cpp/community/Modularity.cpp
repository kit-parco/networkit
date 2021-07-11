// no-networkit-format
/*
 * Modularity.cpp
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt
 */


#include <cmath>
#include <stdexcept>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/Coverage.hpp>
#include <networkit/community/Modularity.hpp>

namespace NetworKit {

Modularity::Modularity() : QualityMeasure(), gTotalEdgeWeight(0.0) {}

void Modularity::setTotalEdgeWeight(double totalEdgeWeight) {
    gTotalEdgeWeight = totalEdgeWeight;
}


double Modularity::getQuality(const Partition& zeta, const Graph& G) {
    assert (G.numberOfNodes() <= zeta.numberOfElements());

    Coverage coverage;
    double cov = coverage.getQuality(zeta, G); // deprecated: intraEdgeWeightSum / gTotalEdgeWeight;
    double expCov; // term $\frac{ \sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }$
    double modularity; // mod = coverage - expected coverage
    if (gTotalEdgeWeight == 0.0) {
        gTotalEdgeWeight = G.totalEdgeWeight(); // compute total edge weight in G
        DEBUG("computed total edge weight: " , gTotalEdgeWeight);
    }

    if (gTotalEdgeWeight == 0.0) {
        ERROR("G: m=" , G.numberOfEdges() , "n=" , G.numberOfNodes());
        throw std::invalid_argument("Modularity is undefined for graphs without edges (including self-loops).");
    }

    std::vector<double> incidentWeightSum(zeta.upperBound(), 0.0); //!< cluster -> sum of the weights of incident edges for all nodes

    // compute volume of each cluster
    G.parallelForNodes([&](node v) {
        // add to cluster weight
        index c = zeta[v];
        assert (zeta.lowerBound() <= c);
        assert (c < zeta.upperBound());
#pragma omp atomic
        incidentWeightSum[c] += G.weightedDegree(v) + G.weight(v,v); // account for self-loops a second time
    });

    // compute sum of squared cluster volumes and divide by squared graph volume
    // double totalIncidentWeight = 0.0; //!< term $\sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 $
    expCov = 0.0;
// double divisor = 4 * totalEdgeWeight * totalEdgeWeight;
// assert (divisor != 0); // do not divide by 0

    #pragma omp parallel for reduction(+:expCov)
    for (omp_index c = static_cast<omp_index>(zeta.lowerBound());
         c < static_cast<omp_index>(zeta.upperBound()); ++c) {
        expCov += ((incidentWeightSum[c] / gTotalEdgeWeight) * (incidentWeightSum[c] / gTotalEdgeWeight )) / 4; // squared
    }

    DEBUG("expected coverage: " , expCov);

    // assert ranges of coverage
    assert(cov <= 1.0);
    assert(cov >= 0.0);
    assert(expCov <= 1.0);
    assert(expCov >= 0.0);

    modularity = cov - expCov;
    DEBUG("modularity = " , modularity);

    // reset totalEdgeWeight
    gTotalEdgeWeight = 0.0;

    assert(! std::isnan(modularity)); // do not return NaN
    // do not return anything not in the range of modularity values
    assert(modularity >= -0.5);
    assert(modularity <= 1);
    return modularity;
}

} /* namespace NetworKit */
