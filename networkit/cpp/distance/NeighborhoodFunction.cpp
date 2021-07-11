// no-networkit-format
/*
* NeighborhoodFunction.cpp
*
*  Created on: 30.03.2016
*      Author: Maximilian Vogel
*/
#include <cmath>
#include <map>
#include <omp.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/NeighborhoodFunction.hpp>
#include <networkit/graph/BFS.hpp>

namespace NetworKit {

NeighborhoodFunction::NeighborhoodFunction(const Graph& G) : Algorithm(), G(&G) {
    if (G.isDirected())
        throw std::runtime_error("current implementation can only deal with undirected graphs");
    ConnectedComponents cc(G);
    cc.run();
    if (cc.numberOfComponents() > 1) 
        throw std::runtime_error("current implementation only runs on graphs with 1 connected component");
}

void NeighborhoodFunction::run() {
    count max_threads = (count)omp_get_max_threads();
    std::vector<std::map<count, count>> nf(max_threads);
    G->parallelForNodes([&](node u){
        index tid = omp_get_thread_num();
        Traversal::BFSfrom(*G, u, [&](node, count dist) {
            nf[tid][dist] += 1;
        });
    });
    count size = 0;
    for (index i = 0; i < max_threads; ++i) {
        size = std::max(size, (count)nf[i].size());
    }
    result = std::vector<count>(size-1, 0);
    for (const auto& local_nf : nf) {
        for (const auto& elem : local_nf) {
            if (elem.first > 0) {
                result[elem.first-1] += elem.second;
            }
        }
    }
    for (index i = 1; i < size-1; ++i) {
        result[i] += result[i-1];
    }
    hasRun = true;
}

std::vector<count> NeighborhoodFunction::getNeighborhoodFunction() const {
    if(!hasRun) {
        throw std::runtime_error("Call run()-function first.");
    }
    return result;
}

} // namespace NetworKit
