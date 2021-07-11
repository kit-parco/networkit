// no-networkit-format
/*
 * Coverage.cpp
 *
 *  Created on: 02.02.2013
 *      Author: Christian Staudt
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/Coverage.hpp>

namespace NetworKit {

double Coverage::getQuality(const Partition& zeta, const Graph& G) {

    double cov = 0.0; // term $\frac{\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)}{\sum_{e \in E} \omega(e)}$
    double totalEdgeWeight = G.totalEdgeWeight(); // add edge weight
    TRACE("total edge weight = " , totalEdgeWeight);

    if (totalEdgeWeight == 0.0) {
        throw std::invalid_argument(
                "Coverage is undefined for graphs without edges (including self-loops).");
    }

    std::vector<double> intraEdgeWeight(zeta.upperBound(), 0.0); // cluster -> weight of its internal edges

    // compute intra-cluster edge weights per cluster
    // TODO: Make parallel, protect intraEdgeWeight[c]
    G.forEdges(
            [&](node u, node v, edgeweight ew) {
                assert (u < zeta.numberOfElements());
                assert (v < zeta.numberOfElements());
                index c = zeta[u];
                index d = zeta[v];
                if (c == d) {
#ifdef DEBUG
                    if ((c >= zeta.upperBound()) || (c < zeta.lowerBound())) {
                        ERROR("c=" , c , " = zeta(" , u , ") is larger than upper bound: " , zeta.upperBound());
                    }
#endif
                    assert ((zeta.lowerBound()) <= c && (c < zeta.upperBound()));
                    intraEdgeWeight[c] += ew;
                } // else ignore edge
            });

    double intraEdgeWeightSum = 0.0; //!< term $\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)$
    #pragma omp parallel for reduction(+:intraEdgeWeightSum)
    for (omp_index c = static_cast<omp_index>(zeta.lowerBound());
         c < static_cast<omp_index>(zeta.upperBound()); ++c) {
        intraEdgeWeightSum += intraEdgeWeight[c];
    }
    TRACE("total intra-cluster edge weight = " , intraEdgeWeightSum);

    cov = intraEdgeWeightSum / totalEdgeWeight;
    DEBUG("coverage = " , cov);

    assert(cov <= 1.0);
    assert(cov >= 0.0);
    return cov;
}

} /* namespace NetworKit */
