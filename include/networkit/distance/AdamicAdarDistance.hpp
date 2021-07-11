// no-networkit-format
/*
 * AdamicAdarDistance.hpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#ifndef NETWORKIT_DISTANCE_ADAMIC_ADAR_DISTANCE_HPP_
#define NETWORKIT_DISTANCE_ADAMIC_ADAR_DISTANCE_HPP_

#include <networkit/distance/NodeDistance.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * An implementation of the Adamic Adar distance measure.
 */
class AdamicAdarDistance final : public NodeDistance {

private:
    std::vector<double> aaDistance; //result vector

public:

    /**
     * @param G The graph.
     */
    AdamicAdarDistance(const Graph& G);

    /**
     * Computes the Adamic Adar distances of all connected pairs of nodes.
     * REQ: Needs to be called before distance() and getEdgeScores() deliver meaningful results!
     */
     void preprocess() override;

    /**
     * Returns the Adamic Adar distance between node @a u and node @a v.
     * @return Adamic Adar distance between the two nodes.
     */
     double distance(node u, node v) override;

    /**
     * Returns the Adamic Adar distances between all connected nodes.
     * @return Vector containing the Adamic Adar distances between all connected pairs of nodes.
     */
     std::vector<double> getEdgeScores() override;

};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_ADAMIC_ADAR_DISTANCE_HPP_
