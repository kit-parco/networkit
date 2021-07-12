/*
 * ParallelPartitionCoarsening.cpp
 *
 *  Created on: 28.01.2014
 *      Author: cls
 */

// networkit-format

#include <numeric>
#include <omp.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/graph/GraphBuilder.hpp>

namespace NetworKit {

ParallelPartitionCoarsening::ParallelPartitionCoarsening(const Graph &G, const Partition &zeta,
                                                         bool parallel)
    : GraphCoarsening(G), zeta(zeta), parallel(parallel) {}

void ParallelPartitionCoarsening::run() {
    Partition nodeToSuperNode = zeta;
    nodeToSuperNode.compact((zeta.upperBound() <= G->upperNodeIdBound()));
    index numParts = nodeToSuperNode.upperBound();

    // sort fine vertices by coarse vertices
    std::vector<index> partBegin(numParts + 2, 0);
    std::vector<node> nodesSortedByPart(G->numberOfNodes());
    G->forNodes([&](const node u) { partBegin[nodeToSuperNode[u] + 2]++; });
    std::partial_sum(partBegin.begin(), partBegin.end(), partBegin.begin());
    G->forNodes([&](const node u) { nodesSortedByPart[partBegin[nodeToSuperNode[u] + 1]++] = u; });

    Gcoarsened = Graph(numParts, true, false);

    auto aggregateEdgeWeights = [&](node su, count &numEdges, count &numSelfLoops,
                                    std::vector<edgeweight> &incidentPartWeights,
                                    std::vector<node> &incidentParts) {
        for (index i = partBegin[su]; i < partBegin[su + 1]; ++i) {
            node u = nodesSortedByPart[i];
            G->forNeighborsOf(u, [&](node v, edgeweight ew) {
                const node sv = nodeToSuperNode[v];
                if (sv != su || u >= v) {
                    if (incidentPartWeights[sv] == 0.0) {
                        incidentParts.push_back(sv);
                    }
                    incidentPartWeights[sv] += ew;
                }
            });
        }

        numEdges += incidentParts.size();
        if (incidentPartWeights[su] != 0.0) {
            numSelfLoops += 1;
            numEdges -= 1;
        }

        Gcoarsened.preallocateUndirected(su, incidentParts.size());
        for (node sv : incidentParts) {
            Gcoarsened.addPartialEdge(unsafe, su, sv, incidentPartWeights[sv]);
            incidentPartWeights[sv] = 0.0;
        }
        incidentParts.clear();
    };

    count numEdges = 0;
    count numSelfLoops = 0;
    if (!parallel) {
        // The code has duplication because even the overhead of opening the parallel section was
        // too much if used on lots of very small graphs, as in Ego-Splitting for example
        std::vector<edgeweight> incidentPartWeights(numParts, 0.0);
        std::vector<node> incidentParts;
        incidentParts.reserve(numParts);
        for (node su = 0; su < numParts; ++su) {
            aggregateEdgeWeights(su, numEdges, numSelfLoops, incidentPartWeights, incidentParts);
        }
    } else {
#pragma omp parallel
        {
            std::vector<edgeweight> incidentPartWeights(numParts, 0.0);
            std::vector<node> incidentParts;
            incidentParts.reserve(numParts);

            count localEdges = 0;
            count localSelfLoops = 0;

#pragma omp for schedule(guided) nowait
            for (omp_index su = 0; su < static_cast<omp_index>(numParts); ++su) {
                aggregateEdgeWeights(su, localEdges, localSelfLoops, incidentPartWeights,
                                     incidentParts);
            }

#pragma omp atomic
            numEdges += localEdges;

#pragma omp atomic
            numSelfLoops += localSelfLoops;
        }
    }

    Gcoarsened.setNumberOfSelfLoops(unsafe, numSelfLoops);
    Gcoarsened.setEdgeCount(unsafe, numEdges / 2 + numSelfLoops);

    this->nodeMapping = nodeToSuperNode.getVector();
    hasRun = true;
}

} /* namespace NetworKit */
