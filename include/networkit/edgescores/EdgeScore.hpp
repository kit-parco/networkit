/*
 * EdgeScore.hpp
 *
 *  Created on: 18.08.2015
 *      Author: Gerd Lindner
 */

#ifndef NETWORKIT_EDGESCORES_EDGE_SCORE_HPP_
#define NETWORKIT_EDGESCORES_EDGE_SCORE_HPP_

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
/**
 * Abstract base class for an edge score.
 */
template <typename T>
class EdgeScore : public Algorithm {

public:
    EdgeScore(const Graph &G);

    /** Compute the edge score. */
    void run() override;

    /** Get a vector containing the score for each edge in the graph.
     *
     * @return the edge scores calculated by @ref run().
     */
    virtual std::vector<T> scores() const;

    /** Get the edge score of the edge with the given edge id.
     */
    virtual T score(edgeid eid);

    /** Get the edge score of the given edge.
     */
    virtual T score(node u, node v);

protected:
    const Graph *G;
    std::vector<T> scoreData;
};

} // namespace NetworKit

#endif // NETWORKIT_EDGESCORES_EDGE_SCORE_HPP_
