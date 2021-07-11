/*
 * TopHarmonicCloseness.hpp
 *
 * Created on: 25.02.2018
 *		Authors: nemes, Eugenio Angriman
 */

#ifndef NETWORKIT_CENTRALITY_TOP_HARMONIC_CLOSENESS_HPP_
#define NETWORKIT_CENTRALITY_TOP_HARMONIC_CLOSENESS_HPP_

#include <memory>

#include <networkit/base/Algorithm.hpp>
#include <networkit/components/WeaklyConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class TopHarmonicCloseness final : public Algorithm {
public:
    /**
     * Finds the top k nodes with highest harmonic closeness centrality faster
     * than computing it for all nodes. The implementation is based on
     * "Computing Top-k Centrality Faster in Unweighted Graphs", Bergamini et
     * al., ALENEX16. The algorithm also works with weighted graphs but only
     * with the NBcut variation. We recommend to use useNBbound = false for
     * (weighted) complex networks (or networks with small diameter) and
     * useNBbound = true for street networks (or networks with large
     * diameters). Notice that the worst case running time of the algorithm is
     * O(nm), where n is the number of nodes and m is the number of edges.
     * However, for most real-world networks the empirical running time is
     * O(m).
     *
     * @param G The input graph. If useNBbound is set to true, edge weights will be ignored.
     * @param k Number of nodes with highest harmonic closeness that have to be found.
     * For example, k = 10, the top 10 nodes with highest harmonic closeness will be computed.
     * @param useNBbound If true, the NBbound variation will be used, otherwise NBcut.
     */
    explicit TopHarmonicCloseness(const Graph &G, count k = 1, bool useNBbound = false);

    ~TopHarmonicCloseness() override;

    /**
     * Computes top-k harmonic closeness on the graph passed in the constructor.
     */
    void run() override;

    /**
     * Returns a list with the k nodes with highest harmonic closeness.
     * WARNING: closeness centrality of some nodes below the top-k could be equal
     * to the k-th closeness, we call them trail. Set the parameter includeTrail
     * to true to also include those nodes but consider that the resulting vector
     * could be longer than k.
     *
     * @param includeTrail Whether or not to include trail nodes.
     * @return The list of the top-k nodes.
     */
    std::vector<node> topkNodesList(bool includeTrail = false);

    /**
     * Returns a list with the scores of the k nodes with highest harmonic
     * closeness.
     * WARNING: closeness centrality of some nodes below the top-k could
     * be equal to the k-th closeness, we call them trail. Set the parameter
     * includeTrail to true to also include those centrality values but consider
     * that the resulting vector could be longer than k.
     *
     * @param includeTrail Whether or not to include trail centrality value.
     * @return The closeness centralities of the k most central nodes.
     */
    std::vector<edgeweight> topkScoresList(bool includeTrail = false);

private:
    const Graph *G;
    const count k;
    const bool useNBbound;

    std::vector<double> hCloseness;
    std::vector<count> reachableNodes;
    std::vector<std::vector<uint8_t>> visitedGlobal;
    std::vector<uint8_t> tsGlobal;

    std::vector<node> topKNodes;
    std::vector<double> topKScores;

    // For NBbound
    std::vector<count> reachU, reachL;
    std::vector<std::vector<count>> numberOfNodesAtLevelGlobal;
    std::vector<std::vector<node>> nodesAtCurrentLevelGlobal;
    std::vector<std::vector<std::vector<node>>> nodesAtLevelGlobal;
    std::vector<double> levelImprovement;
    std::unique_ptr<WeaklyConnectedComponents> wccPtr;
    void updateTimestamp() {
        auto &visited = visitedGlobal[omp_get_thread_num()];
        auto &ts = tsGlobal[omp_get_thread_num()];
        if (ts++ == std::numeric_limits<uint8_t>::max()) {
            ts = 1;
            std::fill(visited.begin(), visited.end(), 0);
        }
    }

    template <class Type>
    struct Greater {
        Greater(const std::vector<Type> &vec) : vec(&vec) {}
        bool operator()(node u, node v) const noexcept { return (*vec)[u] > (*vec)[v]; }
        const std::vector<Type> *vec;
    };

    template <class Type>
    struct Less {
        Less(const std::vector<Type> &vec) : vec(&vec) {}
        bool operator()(node u, node v) const noexcept { return (*vec)[u] < (*vec)[v]; }
        const std::vector<Type> *vec;
    };

    tlx::d_ary_addressable_int_heap<node, 2, Greater<double>> prioQ;
    tlx::d_ary_addressable_int_heap<node, 2, Less<double>> topKNodesPQ;
    std::vector<node> trail;

    // For weighted graphs
    std::vector<tlx::d_ary_addressable_int_heap<node, 2, Less<edgeweight>>> dijkstraHeaps;
    std::vector<std::vector<edgeweight>> distanceGlobal;
    edgeweight minEdgeWeight;

    omp_lock_t lock;

    void runNBcut();
    void runNBbound();
    bool bfscutUnweighted(node source, double kthCloseness);
    bool bfscutWeighted(node source, double kthCloseness);
    void bfsbound(node source);
    void computeReachableNodes();
    void computeReachableNodesBounds();
    void computeNeighborhoodBasedBound();
    void updateTopkPQ(node u);

    double initialBoundNBcutUnweighted(node u) const noexcept;
    double initialBoundNBcutWeighted(node u) const noexcept;

    template <class Type>
    std::vector<Type> resizedCopy(const std::vector<Type> &vec) const noexcept;
};

} // namespace NetworKit
#endif // NETWORKIT_CENTRALITY_TOP_HARMONIC_CLOSENESS_HPP_
