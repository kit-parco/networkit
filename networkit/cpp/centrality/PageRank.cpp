/*
 * PageRank.cpp
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/centrality/PageRank.hpp>

namespace NetworKit {

PageRank::PageRank(const Graph &G, double damp, double tol, bool normalized)
    : Centrality(G, true), damp(damp), tol(tol), normalized(normalized) {}

void PageRank::run() {
    Aux::SignalHandler handler;
    const auto n = G.numberOfNodes();
    const auto z = G.upperNodeIdBound();

    const auto teleportProb = (1.0 - damp) / static_cast<double>(n);
    scoreData.resize(z, 1.0 / static_cast<double>(n));
    std::vector<double> pr = scoreData;

    std::vector<double> deg(z, 0.0);
    G.parallelForNodes([&](const node u) { deg[u] = static_cast<double>(G.weightedDegree(u)); });

    iterations = 0;

    auto sumL1Norm = [&](const node u) { return std::abs(scoreData[u] - pr[u]); };

    auto sumL2Norm = [&](const node u) {
        const auto d = scoreData[u] - pr[u];
        return d * d;
    };

    auto converged([&]() {
        if (iterations >= maxIterations) {
            return true;
        }

        if (norm == Norm::L2Norm) {
            return std::sqrt(G.parallelSumForNodes(sumL2Norm)) <= tol;
        }

        return G.parallelSumForNodes(sumL1Norm) <= tol;
    });

    bool isConverged = false;
    do {
        handler.assureRunning();
        G.balancedParallelForNodes([&](const node u) {
            pr[u] = 0.0;
            G.forInEdgesOf(u, [&](const node u, const node v, const edgeweight w) {
                // note: inconsistency in definition in Newman's book (Ch. 7) regarding directed
                // graphs we follow the verbal description, which requires to sum over the incoming
                // edges
                pr[u] += scoreData[v] * w / deg[v];
            });
            pr[u] *= damp;
            pr[u] += teleportProb;
        });

        ++iterations;
        isConverged = converged();
        std::swap(pr, scoreData);
    } while (!isConverged);

    handler.assureRunning();

    if (!normalized) {
        const auto sum = G.parallelSumForNodes([&](const node u) { return scoreData[u]; });

        // make sure scoreData sums up to 1
        assert(!Aux::NumericTools::equal(sum, 0.0, 1e-15));
        G.parallelForNodes([&](const node u) { scoreData[u] /= sum; });
    } else {
        // calculate sum of the Pagerank of dangling Nodes for normalization
        const auto sum = G.parallelSumForNodes([&](const node u) {
            if (G.degree(u) == 0.0)
                return scoreData[u];
            return 0.0;
        });

        const auto normFactor = (1.0 / n) * ((1 - damp) + (damp * sum));

        G.parallelForNodes([&](const node u) { scoreData[u] /= normFactor; });
    }

    // calculate the maxium
    if (!normalized) {
        max = 1.0; // upper bound, could be tighter by assuming e.g. a star graph with n nodes
    } else {
        max = scoreData[0]; // Unlike regular Page Rank there isn't a universal upper bound for
                            // normalized page rank. So we need to iterate.
        G.balancedParallelForNodes(
            [&](const node u) { Aux::Parallel::atomic_max(max, scoreData[u]); });
    }
    hasRun = true;
}

double PageRank::maximum() {
    return max.load(std::memory_order_relaxed);
}

} /* namespace NetworKit */
