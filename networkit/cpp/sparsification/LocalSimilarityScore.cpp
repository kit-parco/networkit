// no-networkit-format
/*
 * LocalSimilarityScore.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#include <cmath> //log
#include <set>
#include <networkit/sparsification/LocalSimilarityScore.hpp>

namespace NetworKit {

LocalSimilarityScore::LocalSimilarityScore(const Graph& G, const std::vector<count>& triangles) :
    EdgeScore<double>(G), triangles(&triangles) {}

void LocalSimilarityScore::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    /*
     * For each edge, we calculate the minimum required sparsification exponent e
     * such that the edge is contained in the sparse graph.
     */

    std::vector<double> sparsificationExp(G->upperEdgeIdBound(), 0.0);

    // G->parallelForNodes is replaced here by it's code-logic due to bugs in clang 12.0.0 and 12.0.1
    #pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(G->upperNodeIdBound()); ++i) {
        if(!G->hasNode(i)) {
            continue;
        }
        count d = G->degree(i);

        /* The top d^e edges (sorted by similarity)
         * are to be kept in the graph. */

        std::vector<AttributizedEdge<double>> neighbors;
        neighbors.reserve(G->degree(i));
        G->forNeighborsOf(i, [&](node, node j, edgeid eid) {
            double sim = static_cast<double>((*triangles)[eid]) * 1.0 / static_cast<double>(G->degree(i) + G->degree(j) - (*triangles)[eid]);
            neighbors.emplace_back(i, j, eid, sim);
        });
        std::sort(neighbors.begin(), neighbors.end());

        count rank = 1;

        /**
         * By convention, we want to the edges with highest "similarity" or "cohesion" to have values close to 1,
         * so we invert the range.
         */

        #pragma omp critical
        for(std::vector<AttributizedEdge<double>>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
            edgeid eid = it->eid;

            double e = 1.0; //If the node has only one neighbor, the edge will be kept anyway.
            if (d > 1)
                e = 1 - (log(rank) / log(d));

            sparsificationExp[eid] = std::max(e, sparsificationExp[eid]);
            rank++;
        }

    };

    scoreData = std::move(sparsificationExp);
    hasRun = true;
}

double LocalSimilarityScore::score(node, node) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

double LocalSimilarityScore::score(edgeid) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
