// no-networkit-format
/*
 * STSP.hpp
 *
 *  Created on: 15.06.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_DISTANCE_STSP_HPP_
#define NETWORKIT_DISTANCE_STSP_HPP_

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Abstract base class for source-target shortest path algorithms.
 */
class STSP : public Algorithm {

public:
    /**
     * Creates the STSP class for a graph @a G, source node @a source, and
     * target node @a target.
     *
     * @param G The graph.
     * @param source The source node.
     * @param target The target node.
     * @param storePred If true, the algorithm will also store the predecessors
     * and reconstruct a shortest path from @a source and @a target.
     */
    STSP(const Graph &G, node source, node target, bool storePred = true)
        : G(&G), source(source), target(target), storePred(storePred) {
        if (!G.hasNode(source))
            throw std::runtime_error("Error: source node not in the graph!");
        if (!G.hasNode(target))
            throw std::runtime_error("Error: target node not in the graph!");
        if (source == target)
            INFO("Source and target nodes are equal!");
    }

    /**
     * Returns a shortest path from the source node to the target node (without
     * including them). Note: the shortest path can be constructed only if the
     * algorithm is executed with @a storePred set to true.
     *
     * @return A shortest path from the @a source to @a target.
     */
    virtual std::vector<node> getPath() const {
        checkStorePredecessors();
        assureFinished();
        return path;
    }

    /**
     * Returns the predecessor nodes from the @a target to the source. Note:
     * predecessors are stored only if the algorithm is executed with @a
     * storePred set to true.
     *
     * @return The list of predecessors from @a target to @a source.
     */
    const std::vector<node>& getPredecessors() const {
        checkStorePredecessors();
        assureFinished();
        return pred;
    }

    /**
     * Returns the distance from the source node to the target node
     * @return The distance from source to the target node.
     */
    virtual edgeweight getDistance() const = 0;

protected:
    const Graph *G;
    node source, target;
    const bool storePred;
    std::vector<node> pred, path;

    // Builds a source-target shortest path using the predecessors
    void buildPath();

    void checkStorePredecessors() const {
        if (!storePred)
            WARN("Predecessors not stored.");
    }

    void init() {
        hasRun = false;
        if (storePred)
            pred.resize(G->upperNodeIdBound());
    }
};

inline void STSP::buildPath() {
    path.clear();
    node t = target;
    while (pred[t] != source) {
        path.push_back(pred[t]);
        t = pred[t];
    }
    std::reverse(path.begin(), path.end());
}

} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_STSP_HPP_
