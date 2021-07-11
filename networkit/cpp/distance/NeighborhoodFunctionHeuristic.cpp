// no-networkit-format
/*
* NeighborhoodFunctionHeuristic.cpp
*
*      Author: Maximilian Vogel
*/

#include <cmath>
#include <map>
#include <omp.h>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/distance/NeighborhoodFunctionHeuristic.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

NeighborhoodFunctionHeuristic::NeighborhoodFunctionHeuristic(const Graph& G, count nSamples, SelectionStrategy strategy) :
    Algorithm(),
    G(&G),
    nSamples(!nSamples ? (count)ceil(std::max((double)0.15f * G.numberOfNodes(), sqrt(G.numberOfEdges()))) : nSamples),
    strategy(strategy) {

    if (G.isDirected())
        throw std::runtime_error("current implementation can only deal with undirected graphs");
    ConnectedComponents cc(G);
    cc.run();
    if (cc.numberOfComponents() > 1)
        throw std::runtime_error("current implementation only runs on graphs with 1 connected component");

    if (strategy != SPLIT && strategy != RANDOM) {
        throw std::runtime_error("unknown strategy, choose either split or random");
    }
}

void NeighborhoodFunctionHeuristic::run() {
    count maxThreads = (count)omp_get_max_threads();
    count dia;
    if (!G->isWeighted()) {
        Diameter diam(*G);
        diam.run();
        dia = diam.getDiameter().first;
    } else {
        Graph Gcopy = GraphTools::toUnweighted(*G);
        Diameter diam(Gcopy);
        diam.run();
        dia = diam.getDiameter().first;
    }

    std::vector<node> start_nodes(nSamples);
    if (strategy == SPLIT) {
        start_nodes = split(*G, nSamples);
    } else if (strategy == RANDOM) {
        start_nodes = random(*G, nSamples);
    }

    // run nSamples BFS and count the distances.
    std::vector<std::vector<count>> nf(maxThreads, std::vector<count>(dia+1,0));
    #pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(nSamples); ++i) {
        count tid = omp_get_thread_num();
        node u = start_nodes[i];
        Traversal::BFSfrom(*G, u, [&](node, count dist) {
            nf[tid][dist] += 1;
        });
    }

    count n = G->numberOfNodes();
    result.resize(dia, 0);
    count start = 1;
    count end = dia;
    // enhancements to the result
    if (true) {
        count m = G->numberOfEdges();
        result[0] = 2 * m;
        result.back() = n * (n-1);
        start += 1;
        end -= 1;
    }
    // accumulate thread local results and the nf
    double norm_factor = n / (double)nSamples;
    for (index dist = start; dist <= end; ++dist) {
        // accumulate thread local results for each distance
        count tmp = 0;
        #pragma omp parallel for reduction(+:tmp)
        for (omp_index tid = 0; tid < static_cast<omp_index>(nf.size()); ++tid) {
            tmp += nf[tid][dist];
        }
        // accumulate nf
        result[dist-1] = round(tmp * norm_factor);
        if (dist > 1) {
            result[dist-1] += result[dist-2];
        }
    }

    hasRun = true;
}

std::vector<count> NeighborhoodFunctionHeuristic::getNeighborhoodFunction() const {
    if(!hasRun) {
        throw std::runtime_error("Call run()-function first.");
    }
    return result;
}

std::vector<node> NeighborhoodFunctionHeuristic::random(const Graph& G, count nSamples) {
    std::vector<node> start_nodes(nSamples, 0);
    // the vector of start nodes is chosen completely at random with the graphs "randomNode()" function.
    for (index i = 0; i < nSamples; ++i) {
        start_nodes[i] = GraphTools::randomNode(G);
    }
    return start_nodes;
}

std::vector<node> NeighborhoodFunctionHeuristic::split(const Graph& G, count nSamples) {
    // sort the nodes by node degree or sum of degrees of its neighbors
    std::vector<count> nodeDeg(G.upperNodeIdBound(), none);
    G.parallelForNodes([&](node u) {
            nodeDeg[u] = G.degree(u);
    });
    std::vector<node> nodes;
    nodes.reserve(G.numberOfNodes());
    G.forNodes([&](node u) { nodes.push_back(u); });
    nodes.erase(std::remove_if(nodes.begin(), nodes.end(), [](node u){return u == none;}), nodes.end());
    //std::random_shuffle(nodes.begin(), nodes.end());
    Aux::Parallel::sort(nodes.begin(), nodes.end(), [&nodeDeg](const node& a, const node& b) {return nodeDeg[a] < nodeDeg[b];});
    std::vector<node> start_nodes(nSamples, 0);
    // every (n/nSamples)-th node is selected as start node
    auto stepwidth = G.numberOfNodes()/nSamples;
    for (index i = 0; i < nSamples; ++i) {
        start_nodes[i] = nodes[i * stepwidth];
    }
    return start_nodes;
}

} /* namespace NetworKit */
