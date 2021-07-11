/*
 *  AllSimplePaths.hpp
 *
 *  Created on: 23.06.2017
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_REACHABILITY_ALL_SIMPLE_PATHS_HPP_
#define NETWORKIT_REACHABILITY_ALL_SIMPLE_PATHS_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Determines all the possible simple paths from a given source node to a target node of a directed
 * unweighted graph. It also accepts a cutoff value i.e. the maximum length of paths.
 */
class AllSimplePaths final : public Algorithm {

public:
    /**
     * Creates the AllSimplePaths class for @a G, source @a s and target @a t.
     *
     * @param G The graph.
     * @param source The source node.
     * @param target The target node.
     * @param cutoff The maximum length of the paths.
     */
    AllSimplePaths(const Graph &G, node source, node target, count cutoff = none);

    ~AllSimplePaths() override = default;

    /**
     * This method computes all possible paths from a given source node to a target node.
     */
    void run() override;

    /**
     * This method returns the number of simple paths from the source node to the target node.
     */
    count numberOfSimplePaths();

    /*
     * This method returns a vector that contains all the simple paths from a source node to a
     * target node represented by vectors. Each path contains the source node as the first element
     * and the target node as the last element.
     */
    std::vector<std::vector<node>> getAllSimplePaths();

    /*
     * This method iterates over all the simple paths and it is far more efficient than calling
     * getAllSimplePaths().
     */
    template <typename L>
    void forAllSimplePaths(L handle);

    /*
     * This method iterates in parallel over all the simple paths and it is far more efficient than
     * calling getAllSimplePaths().
     */
    template <typename L>
    void parallelForAllSimplePaths(L handle);

private:
    // This method computes all the paths after a reverse BFS from the target node and a normal BFS
    // from the source node.
    void computePaths();

    // This method returns a queue that contains all the nodes that could be part of a path from the
    // source to the target that crosses @s.
    std::vector<node> getAvailableSources(node s, count pathLength = 0);

    // The graph
    const Graph *G;
    // The source node
    node source;
    // The target node
    node target;
    // The cutoff i.e. maximum length of paths from source to target. It is optional.
    count cutoff;

    // This vector contains the distance from each node to the target node.
    std::vector<count> distanceToTarget;
    // This vector contains the distance from the source node to each node.
    std::vector<count> distanceFromSource;
    // This vector contains all the possible paths from source to target.
    std::vector<std::vector<node>> paths;
};

inline count AllSimplePaths::numberOfSimplePaths() {
    assureFinished();
    return paths.size();
}

inline std::vector<std::vector<node>> AllSimplePaths::getAllSimplePaths() {
    assureFinished();
    return paths;
}

template <typename L>
void AllSimplePaths::forAllSimplePaths(L handle) {
    assureFinished();
    for (std::vector<std::vector<node>>::iterator it = paths.begin(); it != paths.end(); ++it) {
        handle(*it);
    }
}

template <typename L>
void AllSimplePaths::parallelForAllSimplePaths(L handle) {
    assureFinished();
#pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(paths.size()); ++i) {
        handle(paths[i]);
    }
}

} // namespace NetworKit

#endif // NETWORKIT_REACHABILITY_ALL_SIMPLE_PATHS_HPP_
