// no-networkit-format
/*
 * EigenvectorCentrality.cpp
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

#include <cmath>

#include <networkit/centrality/EigenvectorCentrality.hpp>
#include <networkit/auxiliary/NumericTools.hpp>

namespace NetworKit {

EigenvectorCentrality::EigenvectorCentrality(const Graph& G, double tol):
        Centrality(G, true), tol(tol)
{

}

void EigenvectorCentrality::run() {
    std::vector<double> values(G.upperNodeIdBound(), 1.0);
    scoreData = values;

    double length = 0.0;
    double oldLength = 0.0;

    auto converged([tol=tol](double val, double other) -> bool {
        // compute residual
        return (Aux::NumericTools::equal(val, other, tol));
    });

    do {
        oldLength = length;

        // iterate matrix-vector product
        G.parallelForNodes([&](node u) {
            values[u] = 0.0;
            G.forInEdgesOf(u, [&](node u, node v, edgeweight ew) {
                values[u] += ew * scoreData[v];
            });
        });

        // normalize values
        length = G.parallelSumForNodes([&values](node u) {
            return (values[u] * values[u]);
        });
        length = std::sqrt(length);

        assert(! Aux::NumericTools::equal(length, 1e-16));
        G.parallelForNodes([&values, length](node u) {
            values[u] /= length;
        });

        std::swap(scoreData, values);
    } while (! converged(length, oldLength));

    // check sign and correct if necessary
    if (scoreData[0] < 0) {
        G.parallelForNodes([&](node u) {
            scoreData[u] = std::fabs(scoreData[u]);
        });
    }

    hasRun = true;
}

} /* namespace NetworKit */
