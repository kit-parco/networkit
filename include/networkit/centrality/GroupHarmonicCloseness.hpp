/*
 * GroupHarmonicCloseness.hpp
 *
 * Created on: 15.12.2020
 *     Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_HPP_
#define NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_HPP_

#include <memory>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class GroupHarmonicCloseness final : public Algorithm {

public:
    /**
     * Approximation algorithm for the group-harmonic maximization problem. The computed solutions
     * have a guaranteed $\\lambda(1 - \\frac{1}{2e})$ (directed graphs) and
     * $\\lambda(1 - \\frac{1}/{e})/2$ (undirected graphs) approximation ratio,
     * where $\\lambda$ is the ratio
     * between the minimal and the maximal edge weight. The algorithm is the one proposed in
     * Angriman et al., ALENEX 2021. The worst-case running time of this approach is quadratic, but
     * usually much faster in practice.
     *
     * @param G The input graph.
     * @param k Size of the group of nodes.
     */
    GroupHarmonicCloseness(const Graph &G, count k = 1);

    /**
     * Runs the algorithm.
     */
    void run() override;

    /**
     * Returns the computed group.
     *
     * @return The computed group.
     */
    const std::vector<node> &groupMaxHarmonicCloseness() const;

    /**
     * Computes the group-harmonic score of the group of nodes in the given range.
     *
     * @param graph The input Graph.
     * @param first,last The range that contains the vertices in the group.
     *
     * @returns The score of the group of nodes in the given range.
     */
    template <class InputIt>
    static double scoreOfGroup(const Graph &graph, InputIt first, InputIt last);

    class GroupHarmonicClosenessInterface : public Algorithm {
    public:
        std::vector<node> group;
    };

private:
    const bool weighted;
    // Always one between GroupHarmonicClosenessUnweighted and GroupHarmonicClosenessWeighted, see
    // implementation.
    std::unique_ptr<GroupHarmonicClosenessInterface> impl;

    static double scoreOfGroup(const Graph &graph, const std::vector<node> &inputGroup);
};

template <class InputIt>
double GroupHarmonicCloseness::scoreOfGroup(const Graph &graph, InputIt first, InputIt last) {
    return scoreOfGroup(graph, {first, last});
}

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GROUP_HARMONIC_CLOSENESS_HPP_
