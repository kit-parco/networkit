// no-networkit-format
/*
 *  EstimateBetweenness.cpp
 *
 *  Created on: 13.06.2014
 *      Author: Christian Staudt, Elisabetta Bergamini
 */


#include <networkit/centrality/EstimateBetweenness.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/SSSP.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <memory>
#include <omp.h>

namespace NetworKit {

EstimateBetweenness::EstimateBetweenness(const Graph& G, count nSamples, bool normalized, bool parallel_flag) : Centrality(G, normalized), nSamples(nSamples), parallel_flag(parallel_flag) {
}

void EstimateBetweenness::run() {
    hasRun = false;

    Aux::SignalHandler handler;

    std::vector<node> sampledNodes;

    // sample nodes
    for (count i = 0; i <= nSamples; ++i) {
        sampledNodes.push_back(GraphTools::randomNode(G));
    }


    // thread-local scores for efficient parallelism
    count maxThreads = omp_get_max_threads();
    if (!parallel_flag) maxThreads = 1;
    std::vector<std::vector<double> > scorePerThread(maxThreads, std::vector<double>(G.upperNodeIdBound()));


    auto computeDependencies = [&](node s){
        // run single-source shortest path algorithm
        std::unique_ptr<SSSP> sssp;
        if (G.isWeighted()) {
            sssp = std::make_unique<Dijkstra>(G, s, true, true);
        } else {
            sssp = std::make_unique<BFS>(G, s, true, true);
        }
        if (!handler.isRunning()) return;
        sssp->run();
        if (!handler.isRunning()) return;


        // create stack of nodes in non-decreasing order of distance
        auto stack = sssp->getNodesSortedByDistance();

        // compute dependencies and add the contributions to the centrality score
        std::vector<double> dependency(G.upperNodeIdBound(), 0.0);
        for (auto it = stack.rbegin(); it != stack.rend(); ++it) {
            node t = *it;
            if (t == s){
                continue;
            }
            for (node p : sssp->getPredecessors(t)) {
                // TODO: make weighting factor configurable

                // workaround for integer overflow in large graphs
                bigfloat tmp = sssp->numberOfPaths(p) / sssp->numberOfPaths(t);
                double weight;
                tmp.ToDouble(weight);

                dependency[p] += (double(sssp->distance(p)) / sssp->distance(t)) * weight * (1 + dependency[t]);
            }
            scorePerThread[omp_get_thread_num()][t] += dependency[t];
        }
    };


    #pragma omp parallel for if(parallel_flag)
    for (omp_index i = 0; i < static_cast<omp_index>(sampledNodes.size()); ++i) {
        computeDependencies(sampledNodes[i]);
    }

    if (parallel_flag) {
        scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);

        // add up all thread-local values
        for (const auto &local : scorePerThread) {
            G.parallelForNodes([&](node v){
                scoreData[v] += local[v];
            });
        }
    } else {
        scoreData.swap(scorePerThread[0]);
    }


    const count n = G.numberOfNodes();
    const count pairs = (n-2) * (n-1);

    // extrapolate
    G.parallelForNodes([&](node u) {
        scoreData[u] = scoreData[u] * (2 * static_cast<double>(n) / static_cast<double>(nSamples));

        if (normalized) {
            // divide by the number of possible pairs
            scoreData[u] = scoreData[u] / pairs;
        }
    });

    handler.assureRunning();
    hasRun = true;
}


} /* namespace NetworKit */
