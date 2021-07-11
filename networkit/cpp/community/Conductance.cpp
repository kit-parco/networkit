// no-networkit-format
/*
 * Conductance.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include <networkit/community/Conductance.hpp>

namespace NetworKit {

double Conductance::getQuality(const Partition& zeta, const Graph& G) {
    double cond = 0.0;
    double denom = 0.0;

    assert(zeta.numberOfSubsets() == 2);

    if (G.isWeighted()) {
        // compute denominator
        double vol[2] = {0.0, 0.0};

        G.forNodes([&](node v) {
            vol[zeta[v]] += G.weightedDegree(v);
        });

        // check if 2-partition
        assert(Aux::NumericTools::equal(vol[0] + vol[1], 2 * G.totalEdgeWeight()));

        denom = std::min(vol[0], vol[1]);
    }
    else {
        // compute denominator
        count vol[2] = {0, 0};

        G.forNodes([&](node v) {
            vol[zeta[v]] += G.degree(v);
        });

        // check if 2-partition
        assert(vol[0] + vol[1] == 2 * G.totalEdgeWeight());

        denom = (double) std::min(vol[0], vol[1]);
    }

    EdgeCut ec;
    cond = (double) ec.getQuality(zeta, G) / denom;
    return cond;
}

} /* namespace NetworKit */
