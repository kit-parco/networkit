/*
 * DynApproxBetweennssImpl.cpp
 *
 * Author: Charmaine Ndolo <Charmaine.Ndolo@hu-berlin.de>
 */

// networkit-format
#include "DynApproxBetweennessImpl.hpp"
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/centrality/Centrality.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/distance/DynBFS.hpp>
#include <networkit/distance/DynDijkstra.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

DynApproxBetweennessImpl::DynApproxBetweennessImpl(const Graph &G, double epsilon, double delta,
                                                   bool storePredecessors, double universalConstant)
    : Centrality(G, true), epsilon(epsilon), delta(delta), storePreds(storePredecessors),
      universalConstant(universalConstant) {}

count DynApproxBetweennessImpl::getNumberOfSamples() {
    return r;
}

void DynApproxBetweennessImpl::run() {
    if (G.isDirected()) {
        throw std::runtime_error("Invalid argument: G must be undirected.");
    }
    scoreData.clear();
    scoreData.resize(G.upperNodeIdBound());
    u.clear();
    v.clear();
    sampledPaths.clear();

    Diameter diam(G, DiameterAlgo::estimatedPedantic);
    diam.run();
    edgeweight vd = diam.getDiameter().first;

    INFO("estimated diameter: ", vd);
    r = ceil((universalConstant / (epsilon * epsilon)) * (floor(log2(vd - 2)) + 1 - log(delta)));
    INFO("taking ", r, " path samples");
    sssp.clear();
    sssp.resize(r);
    u.resize(r);
    v.resize(r);
    sampledPaths.resize(r);

    for (count i = 0; i < r; i++) {
        DEBUG("sample ", i);
        // sample random node pair
        u[i] = GraphTools::randomNode(G);
        do {
            v[i] = GraphTools::randomNode(G);
        } while (v[i] == u[i]);
        if (G.isWeighted()) {
            sssp[i].reset(new DynDijkstra(G, u[i], storePreds));
        } else {
            sssp[i].reset(new DynBFS(G, u[i], storePreds));
        }
        DEBUG("running shortest path algorithm for node ", u[i]);

        sssp[i]->setTargetNode(v[i]);
        sssp[i]->run();
        if (sssp[i]->distances[v[i]] > 0) { // at least one path between {u, v} exists
            DEBUG("updating estimate for path ", u[i], " <-> ", v[i]);
            // random path sampling and estimation update
            sampledPaths[i].clear();
            node t = v[i];
            while (t != u[i]) {
                // sample z in P_u(t) with probability sigma_uz / sigma_us
                std::vector<std::pair<node, double>> choices;
                if (storePreds) {
                    for (node z : sssp[i]->previous[t]) {
                        // workaround for integer overflow in large graphs
                        bigfloat tmp = sssp[i]->numberOfPaths(z) / sssp[i]->numberOfPaths(t);
                        double weight;
                        tmp.ToDouble(weight);

                        choices.emplace_back(z, weight); // sigma_uz / sigma_us
                    }
                } else {
                    G.forInEdgesOf(t, [&](node t, node z, edgeweight w) {
                        if (Aux::NumericTools::logically_equal(sssp[i]->distances[t],
                                                               sssp[i]->distances[z] + w)) {
                            // workaround for integer overflow in large graphs
                            bigfloat tmp = sssp[i]->numberOfPaths(z) / sssp[i]->numberOfPaths(t);
                            double weight;
                            tmp.ToDouble(weight);

                            choices.emplace_back(z, weight);
                        }
                    });
                }
                DEBUG("Node: ", t);
                DEBUG("Source: ", u[i]);
                assert(!choices.empty());
                node z = Aux::Random::weightedChoice(choices);
                assert(z <= G.upperNodeIdBound());
                if (z != u[i]) {
                    scoreData[z] += 1 / (double)r;
                    sampledPaths[i].push_back(z);
                }
                t = z;
            }
        }
    }

    hasRun = true;
}

void DynApproxBetweennessImpl::update(GraphEvent e) {
    std::vector<GraphEvent> batch(1);
    batch[0] = e;
    updateBatch(batch);
}

void DynApproxBetweennessImpl::updateBatch(const std::vector<GraphEvent> &batch) {
    DEBUG("Updating");
    for (node i = 0; i < r; i++) {
        sssp[i]->updateBatch(batch);
        if (sssp[i]->modified()) {
            // subtract contributions to nodes in the old sampled path
            for (node z : sampledPaths[i]) {
                scoreData[z] -= 1 / (double)r;
            }
            // sample a new shortest path
            sampledPaths[i].clear();
            node t = v[i];
            while (t != u[i]) {
                // sample z in P_u(t) with probability sigma_uz / sigma_us
                std::vector<std::pair<node, double>> choices;
                if (storePreds) {
                    for (node z : sssp[i]->previous[t]) {
                        // workaround for integer overflow in large graphs
                        bigfloat tmp = sssp[i]->numberOfPaths(z) / sssp[i]->numberOfPaths(t);
                        double weight;
                        tmp.ToDouble(weight);

                        choices.emplace_back(z, weight);
                    }
                } else {
                    G.forInEdgesOf(t, [&](node t, node z, edgeweight w) {
                        if (Aux::NumericTools::logically_equal(sssp[i]->distances[t],
                                                               sssp[i]->distances[z] + w)) {
                            // workaround for integer overflow in large graphs
                            bigfloat tmp = sssp[i]->numberOfPaths(z) / sssp[i]->numberOfPaths(t);
                            double weight;
                            tmp.ToDouble(weight);

                            choices.emplace_back(z, weight);
                        }
                    });
                }
                assert(!choices.empty()); // this should fail only if the graph is not connected
                node z = Aux::Random::weightedChoice(choices);
                assert(z <= G.upperNodeIdBound());
                if (z != u[i]) {
                    scoreData[z] += 1 / (double)r;
                    sampledPaths[i].push_back(z);
                }
                t = z;
            }
        }
    }
}
} // namespace NetworKit
