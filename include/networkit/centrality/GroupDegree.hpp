// no-networkit-format
/*
 * GroupDegree.h
 *
 *  Created on: 20.04.2018
 *      Author: Eugenio Angriman
 */

#ifndef NETWORKIT_CENTRALITY_GROUP_DEGREE_HPP_
#define NETWORKIT_CENTRALITY_GROUP_DEGREE_HPP_

#include <omp.h>

#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class GroupDegree : public Algorithm {

public:
    /**
     * Finds the group with the highest group degree centrality according to the
     * definition proposed in 'The centrality of groups and classes' by Everett et
     * al. (The Journal of mathematical sociology, 1999). This is a submodular but
     * non monotone function so the algorithm can find a solution that is at least
     * 1/2 of the optimum. Worst-case running time is quadratic, but usually
     * faster in real-world networks.
     * The 'countGroupNodes' option also count the nodes inside the group in the
     * score, this make the group degree monotone and submodular and the algorithm
     * is guaranteed to return a (1 - 1/e)-approximation of the optimal solution.
     *
     * @param G A graph.
     * @param k Size of the group of nodes
     * @param countGroupNodes if nodes inside the group should be counted in the
     * centrality score.
     */
    GroupDegree(const Graph &G, count k = 1, bool countGroupNodes = true);

    /**
     * Computes the group with maximum degree centrality of the graph passed in
     * the constructor.
     */
    void run() override;

    /**
     * Returns the group with maximum degree centrality.
     */
    std::vector<node> groupMaxDegree();

    /**
     * Returns the score of the group with maximum degree centrality (i.e. the
     * number of nodes outside the group that can be reached in one hop from at
     * least one node in the group).
     */
    count getScore();

    /**
     * Returns the score of the given group.
     */
    count scoreOfGroup(const std::vector<node> &group) const;

protected:
    const Graph &G;
    const count k;
    const bool countGroupNodes;
    count n;
    std::vector<node> group;
    std::vector<int64_t> gain;
    std::vector<bool> reachable;
    std::vector<bool> affected;
    std::vector<bool> inGroup;
    Aux::BucketPQ queue;
    count groupScore;

    void init();
    void updateQueue();
    void updateGroup();
    void computeScore();
    void checkGroup(const std::vector<node> &group) const;
};

inline std::vector<node> GroupDegree::groupMaxDegree() {
    assureFinished();
    return group;
}

inline count GroupDegree::getScore() {
    assureFinished();
    return groupScore;
}

inline void GroupDegree::computeScore() {
    groupScore = std::count(reachable.begin(), reachable.end(), true);

    if (!countGroupNodes) {
        groupScore -= k;
    }
}

inline void GroupDegree::checkGroup(const std::vector<node> &group) const {
    const count z = G.upperNodeIdBound();
    std::vector<bool> check(z, false);
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(group.size()); ++i) {
        node u = group[i];
        if (u >= z) {
            std::stringstream ss;
            ss << "Error: node " << u << " is not in the graph.\n";
            throw std::runtime_error(ss.str());
        }
        if (check[u]) {
            std::stringstream ss;
            ss << "Error: the group contains duplicates of node " << u << ".\n";
            throw std::runtime_error(ss.str());
        }
        check[u] = true;
    }
}

inline count GroupDegree::scoreOfGroup(const std::vector<node> &group) const {
    checkGroup(group);
    std::vector<bool> touched(n, false);
    std::vector<bool> inGroup(n, false);

    for (count i = 0; i < group.size(); ++i) {
        inGroup[group[i]] = true;
    }

    auto processNeighbor = [&](const node u, const node v) {
        if (inGroup[u]) {
            touched[v] = true;
        }
    };

    G.forNodes([&](node v) {
        if (!inGroup[v]) {
            G.forInNeighborsOf(v, [&](node u) { processNeighbor(u, v); });
        }
    });
    count result = std::count(touched.begin(), touched.end(), true);
    if (countGroupNodes) {
        result += group.size();
    }

    return result;
}

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GROUP_DEGREE_HPP_
