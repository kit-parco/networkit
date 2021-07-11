/*
 * TopCloseness.hpp
 *
 *  Created on: 03.10.2014
 *      Author: ebergamini, michele borassi
 */

#ifndef NETWORKIT_CENTRALITY_TOP_CLOSENESS_HPP_
#define NETWORKIT_CENTRALITY_TOP_CLOSENESS_HPP_

#include <memory>

#include <networkit/base/Algorithm.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class TopCloseness final : public Algorithm {
public:
    /**
     * Finds the top k nodes with highest closeness centrality faster than
     * computing it for all nodes, based on "Computing Top-k Closeness Centrality
     * Faster in Unweighted Graphs", Bergamini et al., ALENEX16. The algorithms is
     * based on two independent heuristics, described in the referenced paper. We
     * recommend to use first_heu = true and second_heu = false for complex
     * networks and first_heu = true and second_heu = true for street networks or
     * networks with large diameters. Notice that the worst case running time of
     * the algorithm is O(nm), where n is the number of nodes and m is the number
     * of edges. However, for most networks the empirical running time is O(m).
     *
     * @param G An unweighted graph.
     * @param k Number of nodes with highest closeness that have to be found. For
     * example, if k = 10, the top 10 nodes with highest closeness will be
     * computed.
     * @param first_heu If true, the neighborhood-based lower bound is computed
     * and nodes are sorted according to it. If false, nodes are simply sorted by
     * degree.
     * @param sec_heu If true, the BFSbound is re-computed at each iteration. If
     * false, BFScut is used.
     * @
     */
    TopCloseness(const Graph &G, count k = 1, bool first_heu = true, bool sec_heu = true);

    /**
     * Computes top-k closeness on the graph passed in the constructor.
     */
    void run() override;

    /**
     * Returns a list with the k nodes with highest closeness.
     * WARNING: closeness centrality of some nodes below the top-k could be equal
     * to the k-th closeness, we call them trail. Set the parameter includeTrail
     * to true to also include those nodes but consider that the resulting vector
     * could be longer than k.
     *
     * @param includeTrail Whether or not to include trail nodes.
     */
    std::vector<node> topkNodesList(bool includeTrail = false);

    /**
     * Returns a list with the scores of the k nodes with highest closeness
     * WARNING: closeness centrality of some nodes below the top-k could be equal
     * to the k-th closeness, we call them trail. Set the parameter includeTrail
     * to true to also include those centrality values but consider that the
     * resulting vector could be longer than k.
     *
     * @param includeTrail Whether or not to include trail centrality value.
     */
    std::vector<edgeweight> topkScoresList(bool includeTrail = false);

private:
    const Graph &G;
    count n;
    count k;
    bool first_heu, sec_heu;
    std::vector<node> topk;
    count visEdges;
    count n_op;
    count trail;
    double maxFarness = -1.0;
    count nMaxFarness;
    std::vector<std::vector<count>> nodesPerLevs, sumLevels;
    std::vector<edgeweight> topkScores;
    std::vector<double> farness;
    std::shared_ptr<std::vector<count>> reachLPtr, reachUPtr;

    std::unique_ptr<StronglyConnectedComponents> sccsPtr;

    void init();
    double BFScut(node v, double x, std::vector<bool> &visited, std::vector<count> &distances,
                  std::vector<node> &pred, count &visEdges);
    void computelBound1(std::vector<double> &S);
    void BFSbound(node x, std::vector<double> &S, count &visEdges,
                  const std::vector<bool> &toAnalyze);
    void computeReachable();

    // Returns the node with highest farness
    struct LargerFarness {
        LargerFarness(const std::vector<double> &v_) : v(v_) {}
        bool operator()(node x, node y) const noexcept { return v[x] > v[y]; }

    private:
        const std::vector<double> &v;
    };

    struct SmallerFarness {
        SmallerFarness(const std::vector<double> &v_) : v(v_) {}
        bool operator()(node x, node y) const noexcept { return v[x] < v[y]; }

    private:
        const std::vector<double> &v;
    };
};

inline std::vector<node> TopCloseness::topkNodesList(bool includeTrail) {
    assureFinished();
    if (!includeTrail) {
        auto begin = topk.begin();
        std::vector<node> topkNoTrail(begin, begin + k);
        return topkNoTrail;
    }

    return topk;
}

inline std::vector<edgeweight> TopCloseness::topkScoresList(bool includeTrail) {
    assureFinished();
    if (!includeTrail) {
        auto begin = topkScores.begin();
        std::vector<double> topkScoresNoTrail(begin, begin + k);
        return topkScoresNoTrail;
    }

    return topkScores;
}

} /* namespace NetworKit */
#endif // NETWORKIT_CENTRALITY_TOP_CLOSENESS_HPP_
