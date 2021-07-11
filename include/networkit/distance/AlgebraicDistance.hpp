// no-networkit-format
/*
 * AlgebraicDistance.hpp
 *
 *  Created on: 03.11.2015
 *      Author: Henning Meyerhenke, Christian Staudt, Michael Hamann
 */

#ifndef NETWORKIT_DISTANCE_ALGEBRAIC_DISTANCE_HPP_
#define NETWORKIT_DISTANCE_ALGEBRAIC_DISTANCE_HPP_

#include <networkit/distance/NodeDistance.hpp>


namespace NetworKit {

/**
 * @ingroup distance
 * Algebraic distance assigns a distance value to pairs of nodes
 * according to their structural closeness in the graph.
 * Algebraic distances will become small within dense subgraphs.
 */
class AlgebraicDistance final: public NodeDistance {

public:

    /**
     * @param G The graph.
     * @param numberSystems Number of vectors/systems used for algebraic iteration.
     * @param numberIterations Number of iterations in each system.
     * @param omega attenuation factor influencing convergence speed.
     * @param norm The norm factor of the extended algebraic distance.
     * @param withEdgeScores calculate array of scores for edges {u,v} that equal ad(u,v)
     */
    AlgebraicDistance(const Graph& G, count numberSystems=10, count numberIterations=30, double omega=0.5, index norm=0, bool withEdgeScores=false);

     void preprocess() override;

    /**
     * @return algebraic distance between the two nodes.
     */
     double distance(node u, node v) override;

     std::vector<double> getEdgeScores() override;

private:

    /**
     * initialize vectors randomly
     */
    void randomInit();

    count numSystems; //!< number of vectors/systems used for algebraic iteration
    count numIters; //!< number of iterations in each system
    double omega; //!< attenuation factor influencing the speed of convergence
    index norm;
    const index MAX_NORM = 0;
    bool withEdgeScores;

    std::vector<double> loads; //!< loads[u*numSystems..(u+1)*numSystems]: loads for node u

    std::vector<edgeweight> edgeScores; //!< distance(u,v) for edge {u,v}

};

} /* namespace NetworKit */
#endif // NETWORKIT_DISTANCE_ALGEBRAIC_DISTANCE_HPP_
