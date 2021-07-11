/*
 * SPSP.hpp
 *  Created on: 23.10.2020
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_DISTANCE_SPSP_HPP_
#define NETWORKIT_DISTANCE_SPSP_HPP_

#include <unordered_map>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Class for some-pairs shortest path algorithm.
 */
class SPSP final : public Algorithm {

public:
    /**
     * Creates the SPSP class for @a G.
     *
     * @param G The graph.
     * @param sourcesFirst,sourcesLast Range of the source nodes.
     */
    template <class InputIt>
    explicit SPSP(const Graph &G, InputIt sourcesFirst, InputIt sourcesLast) : G(&G) {
        setSources(sourcesFirst, sourcesLast);
    }

    /**
     * Sets the source nodes.
     *
     * @param sourcesFirst,sourcesLast Range of the new source nodes.
     */
    template <class InputIt>
    void setSources(InputIt sourcesFirst, InputIt sourcesLast) {
        sources.clear();
        sources.insert(sources.begin(), sourcesFirst, sourcesLast);
        sourceIdx.clear();
        for (index i = 0; i < sources.size(); ++i)
            sourceIdx[sources[i]] = i;
    }

    ~SPSP() override = default;

    /**
     * Computes the shortest paths from the source nodes to all other nodes.
     * The algorithm is parallel.
     */
    void run() override;

    /**
     * Returns a vector of weighted distances between all the source nodes to all the other nodes.
     *
     * @return The shortest-path distances from the source nodes to any other node in
     * the graph.
     */
    const std::vector<std::vector<edgeweight>> &getDistances() const {
        assureFinished();
        return distances;
    }

    /**
     * Returns the distance from the source @a u to node @a v or infinity if @a
     * u cannot reach @a v.
     *
     * @param u A source node.
     * @param v A node.
     */
    edgeweight getDistance(node u, node v) const {
        assureFinished();
        return distances[sourceIdx.at(u)][v];
    }

    /**
     * Returns the <node, index> map from source nodes to their index in the vector
     * returned by SPSP::getDistances().
     *
     * @return Map from source nodes to their index in the vector returned by SPSP::getDistances().
     */
    const std::unordered_map<node, index> &getSourceIndexMap() const noexcept { return sourceIdx; }

    /**
     * @return True: the algorithm is parallel.
     */
    bool TLX_DEPRECATED(isParallel() const override) { return true; }

private:
    const Graph *G;
    std::vector<node> sources;
    std::unordered_map<node, index> sourceIdx;
    std::vector<std::vector<edgeweight>> distances;
};

} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_SPSP_HPP_
