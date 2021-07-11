// no-networkit-format
/*
 * KatzCentrality.cpp
 *
 *  Created on: 09.01.2015
 *      Author: Henning
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/centrality/KatzCentrality.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

double KatzCentrality::defaultAlpha(const Graph& G) {
    return 1. / (1. + GraphTools::maxDegree(G));
}

KatzCentrality::KatzCentrality(const Graph& G, double alpha, double beta, double tol):
        Centrality(G, true),
        alpha(alpha == 0 ? KatzCentrality::defaultAlpha(G) : alpha),
        beta(beta), tol(tol)
{
    if (alpha < 0)
        throw std::runtime_error("Alpha must be non negative.");
}

void KatzCentrality::run() {
    std::fill(values.begin(), values.end(), 1.0);
    values.resize(G.upperNodeIdBound(), 1.0);
    scoreData = values;
    double length = 0.0;
    double oldLength = 0.0;

    auto converged = [&](double val, double other) -> bool {
        // compute residual
        return Aux::NumericTools::equal(val, other, tol);
    };

    auto updateScore = [&](node u, node v, edgeweight ew, edgeid) -> void {
        values[u] += ew * alpha * (1 + scoreData[v]);
    };

    do {
        oldLength = length;

        // iterate matrix-vector product
        G.parallelForNodes([&](node u) {
            values[u] = 0.0;

            if (edgeDirection == EdgeDirection::OutEdges)
                G.forNeighborsOf(u, updateScore);
            else if (edgeDirection == EdgeDirection::InEdges)
                G.forInNeighborsOf(u, updateScore);
            else
                throw std::runtime_error("Unsupported edge direction");

            values[u] += beta;
        });

        // normalize values
        length = G.parallelSumForNodes([&](node u) {
            return values[u] * values[u];
        });
        length = std::sqrt(length);

        scoreData = values;
        INFO("oldLength: ", oldLength, ", length: ", length);
    } while (! converged(length, oldLength));

    G.parallelForNodes([&](node u) {
        scoreData[u] /= length;
    });

    hasRun = true;
}

} /* namespace NetworKit */
