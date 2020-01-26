/*
 * ParallelPartitionCoarsening.cpp
 *
 *  Created on: 28.01.2014
 *      Author: cls
 */

#include <numeric>

#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>

namespace NetworKit {

ParallelPartitionCoarsening::ParallelPartitionCoarsening(const Graph &G, const Partition &zeta,
                                                         bool parallel)
        : GraphCoarsening(G), zeta(zeta),
          parallel(parallel) {

}

void ParallelPartitionCoarsening::run() {
    Partition nodeMapping = zeta;
    nodeMapping.compact(nodeMapping.upperBound() <= G->upperNodeIdBound());
    index numParts = nodeMapping.upperBound();

    // Leave out parallel counting sort for now as it requires some more setup.
    std::vector<index> partBegin(numParts + 2, 0);
    std::vector<node> nodesSortedByPart(G->numberOfNodes());
    G->forNodes([&](const node u) {
        partBegin[ nodeMapping[u] + 2 ]++;
    });
    std::partial_sum(partBegin.begin(), partBegin.end(), partBegin.begin());
    G->forNodes([&](const node u) {
        nodesSortedByPart[ partBegin[ nodeMapping[u] + 1 ]++ ] = u;
    });

    Gcoarsened = Graph(numParts, true, false);

    if (!parallel) {

        std::vector<edgeweight> incidentPartWeights(numParts, 0.0);
        std::vector<node> incidentParts;
        incidentParts.reserve(numParts);

        count numEdges = 0;
        count numSelfLoops = 0;
        for (node su = 0; su < numParts; ++su) {
            for (index i = partBegin[su]; i < partBegin[su+1]; ++i) {
                node u = nodesSortedByPart[i];
                G->forNeighborsOf(u, [&](node v, edgeweight ew) {
                    const node sv = nodeMapping[v];
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
        }

        Gcoarsened.m = numEdges / 2 + numSelfLoops;
        Gcoarsened.storedNumberOfSelfLoops = numSelfLoops;

    } else {
        count numEdges = 0;
        count numSelfLoops = 0;

        #pragma omp parallel
        {
            std::vector<edgeweight> incidentPartWeights(numParts, 0.0);
            std::vector<node> incidentParts;
            incidentParts.reserve(numParts);

            count localEdges = 0;
            count localSelfLoops = 0;

            #pragma omp for schedule(guided) nowait
            for (omp_index su = 0; su < static_cast<omp_index>(numParts); ++su) {
                for (index i = partBegin[su]; i < partBegin[su+1]; ++i) {
                    node u = nodesSortedByPart[i];
                    G->forNeighborsOf(u, [&](node v, edgeweight ew) {
                        const node sv = nodeMapping[v];
                        if (sv != su || u >= v) {
                            if (incidentPartWeights[sv] == 0.0) {
                                incidentParts.push_back(sv);
                            }
                            incidentPartWeights[sv] += ew;
                        }
                    });
                }

                localEdges += incidentParts.size();
                if (incidentPartWeights[su] != 0.0) {
                    localSelfLoops += 1;
                    localEdges -= 1;
                }

                Gcoarsened.preallocateUndirected(su, incidentParts.size());

                for (node sv : incidentParts) {
                    Gcoarsened.addPartialEdge(unsafe, su, sv, incidentPartWeights[sv]);
                    incidentPartWeights[sv] = 0.0;
                }

                incidentParts.clear();
            }

            #pragma omp atomic
            numEdges += localEdges;

            #pragma omp atomic
            numSelfLoops += localSelfLoops;
        }

        Gcoarsened.storedNumberOfSelfLoops = numSelfLoops;
        Gcoarsened.m = numEdges / 2 + numSelfLoops;
    }

    this->nodeMapping = nodeMapping.moveVector();
    hasRun = true;
}

} /* namespace NetworKit */
