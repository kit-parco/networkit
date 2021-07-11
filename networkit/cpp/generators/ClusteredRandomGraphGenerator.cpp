/*
 * ClusteredRandomGraphGenerator.cpp
 *
 *  Created on: 28.02.2014
 *      Author: cls
 *              Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NDEBUG
#include <algorithm>
#endif

#include <omp.h>
#include <random>
#include <stdexcept>
#include <unordered_set>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ClusteredRandomGraphGenerator.hpp>
#include <networkit/graph/GraphBuilder.hpp>

namespace NetworKit {

ClusteredRandomGraphGenerator::ClusteredRandomGraphGenerator(count n, count k, double pIntra,
                                                             double pInter)
    : n(n), k(k), pIntra(pIntra), pInter(pInter), zeta(n) {

    if (k == 0)
        throw std::runtime_error("Error: k must be positive");

    const auto isProbability = [](double p) noexcept -> bool { return (p >= 0) && (p <= 1); };

    if (!isProbability(pIntra))
        throw std::runtime_error("Error: pIntra must be in [0, 1]");

    if (!isProbability(pInter))
        throw std::runtime_error("Error: pInter must be in [0, 1]");

    if (pIntra < pInter)
        WARN("Intra-cluster probability is lower than inter-cluster probability.");

    zeta.setUpperBound(k);
}

Graph ClusteredRandomGraphGenerator::generate() {
    GraphBuilder graphBuilder(n);

    // Sequence of all nodes grouped by cluster. Here the vertices of each cluster are stored
    // sequentially one after the other
    std::vector<node> clustersSequence(n);

    // Index in `clustersSequence` where each cluster ends. E.g., clusterEnd[i] is the index in
    // `clustersSequence` of the first vertex in the (i + 1)-th cluster.
    std::vector<size_t> clusterEnd(k);

    {
        // To assign vertices to clusters uniformly at random
        std::uniform_int_distribution<index> clusterDistr{0, k - 1};

        std::vector<count> clusterSizes(k);

        // Assign nodes to clusters
        graphBuilder.parallelForNodes([&zeta = zeta, &clusterDistr, &clusterSizes](node u) {
            const index cluster = clusterDistr(Aux::Random::getURNG());
            zeta.addToSubset(cluster, u);
#pragma omp atomic
            ++clusterSizes[cluster];
        });

        // Last cluster ends at n, each cluster i \in [0, ..., k - 2] ends at
        // n - \sum_{j = i + 1}^{k - 1} cluster[j].size()
        clusterEnd[k - 1] = n;
        for (int64_t i = k - 2; i >= 0; --i)
            clusterEnd[i] = clusterEnd[i + 1] - clusterSizes[i + 1];

        // Counting sort
        graphBuilder.forNodes(
            [&zeta = zeta, &clustersSequence, &clusterEnd, &clusterSizes](node u) {
                const index cluster = zeta.subsetOf(u);
                // Position of u in the sorted array: position where u's cluster ends minus how many
                // vertices are left to be counted in u's cluster.
                const index uIdx = clusterEnd[cluster] - clusterSizes[cluster]--;
                clustersSequence[uIdx] = u;
            });
        assert(std::is_sorted(
            clustersSequence.begin(), clustersSequence.end(),
            [&zeta = zeta](node u, node v) { return zeta.subsetOf(u) < zeta.subsetOf(v); }));
    }

    // Add all the half-out edges from a node u to all the other nodes in the given sequence.
    // To avoid self-loops and multiple edges, only the half-out edges from u to the vertices in
    // (start, end) interval of clustersSequence are added.
    const auto addHalfOutEdges = [&graphBuilder, &clustersSequence](node u, index start, index end,
                                                                    double p) -> void {
        if (end - start <= 1 || p == 0) // No edges to add
            return;

        if (p == 1) { // Add all edges in the interval
            for (index i = start + 1; i < end; ++i)
                graphBuilder.addHalfOutEdge(u, clustersSequence[i]);
            return;
        }

        std::geometric_distribution<index> geom(p);
        do {
            start += geom(Aux::Random::getURNG()) + 1;
            if (start < end)
                graphBuilder.addHalfEdge(u, clustersSequence[start]);
            else
                return;
        } while (true);
    };

#pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
        // Current node
        const node u = clustersSequence[i];
        // Index in `clustersSequence` where the cluster of node `u` ends.
        const index end = clusterEnd[zeta.subsetOf(u)];

        // Add intra-cluster half-out edges.
        addHalfOutEdges(u, i, end, pIntra);

        // Add inter-cluster half-out edges.
        // end - 1 because addHalfOutEdges discards the first vertex in the given sequence, but in
        // this case an edge between u and the first one in the next cluster is possible.
        addHalfOutEdges(u, end - 1, n, pInter);
    }

    return graphBuilder.toGraph(true, true);
}

Partition ClusteredRandomGraphGenerator::getCommunities() {
    return zeta;
}

} // namespace NetworKit
