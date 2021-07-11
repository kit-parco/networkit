// no-networkit-format
/*
 * ApproxGroupBetweenness.hpp
 *
 *  Created on: 13.03.2018
 *      Author: Marvin Pogoda
 */
#ifndef NETWORKIT_CENTRALITY_APPROX_GROUP_BETWEENNESS_HPP_
#define NETWORKIT_CENTRALITY_APPROX_GROUP_BETWEENNESS_HPP_

#include <omp.h>

#include <networkit/base/Algorithm.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class ApproxGroupBetweenness final : public Algorithm {

public:
    /** Constructs the ApproxGroupBetweenness class for a given undirected graph
     * @a G.
     * @param groupSize Size of the set of nodes.
     * @aram epsilon Determines the accuracy of the approximation.
     */
    ApproxGroupBetweenness(const Graph &G, count groupSize, double epsilon);

    /**
     * Approximately computes a set of nodes with maximum groupbetweenness. Based
     * on the algorithm of Mahmoody,Tsourakakis and Upfal.
     */
    void run() override;

    /**
     * Returns a vector of nodes containing the set of nodes with approximated
     * maximum group betweenness.
     */
    std::vector<node> groupMaxBetweenness() const;

    /**
     * Returns the score of the given set.
     */
    double scoreOfGroup(const std::vector<node> &S, bool normalized = false) const;

protected:
    const Graph &G;
    count n;
    std::vector<node> maxGroup;
    const count groupSize;
    const double epsilon;
};

inline std::vector<node> ApproxGroupBetweenness::groupMaxBetweenness() const {
    assureFinished();
    return maxGroup;
}

inline double
ApproxGroupBetweenness::scoreOfGroup(const std::vector<node> &S,
                                     const bool normalized) const {
    if (S.empty())
        throw std::runtime_error("Error: input group is empty");

    count z = G.upperNodeIdBound();
    std::vector<bool> inGroup(z);
    for (node u : S) {
        if (!G.hasNode(u))
            throw std::runtime_error("Error, input group contains nodes not in the graph");
        if (inGroup[u])
            WARN("Input group contains duplicate nodes!");
        inGroup[u] = true;
    }

    std::vector<double> scorePerThread(omp_get_max_threads());
    std::vector<std::vector<double>> deps(omp_get_max_threads(), std::vector<double>(n));
    std::vector<BFS> bfss(omp_get_max_threads(), BFS(G, 0, true, true));

    auto computeDeps = [&](node source) {

        auto &dep = deps[omp_get_thread_num()];
        std::fill(dep.begin(), dep.end(), 0);

        BFS &bfs = bfss[omp_get_thread_num()];
        bfs.setSource(source);
        bfs.run();

        double weight;
        auto sortedNodes = bfs.getNodesSortedByDistance();
        for (auto it = sortedNodes.rbegin(); it != sortedNodes.rend(); ++it) {
            node target = *it;
            for (node pred : bfs.getPredecessors(target)) {
                (bfs.numberOfPaths(pred) / bfs.numberOfPaths(target)).ToDouble(weight);
                if (inGroup[pred]) {
                    if (pred != source) {
                        scorePerThread[omp_get_thread_num()] += dep[pred] + weight * (1 + dep[target]);
                    }
                }
                else {
                    dep[pred] += weight * (1 + dep[target]);
                }
            }
        }
    };

    G.balancedParallelForNodes(computeDeps);

    double result = 0;
    for (double curScore : scorePerThread)
            result += curScore;

    if (normalized) {
        double nPairs = static_cast<double>((z - S.size()) * (z - S.size() - 1));
        if (nPairs <= 0)
            return 0;
        if (!G.isDirected())
            nPairs /= 2;
        result /= nPairs;
    }

    return result;
}

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_APPROX_GROUP_BETWEENNESS_HPP_
