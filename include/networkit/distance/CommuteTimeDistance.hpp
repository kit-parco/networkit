// no-networkit-format
/*
 * CommuteTimeDistance.hpp
 *
 *  Created on: 12.04.2016
 *      Author: ebergamini
 */

#ifndef NETWORKIT_DISTANCE_COMMUTE_TIME_DISTANCE_HPP_
#define NETWORKIT_DISTANCE_COMMUTE_TIME_DISTANCE_HPP_

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/numerics/LAMG/Lamg.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 *
 * CommuteTimeDistance edge centrality.
 *
 */
class CommuteTimeDistance final : public Algorithm {

public:
    /**
     * Constructs the CommuteTimeDistance class for the given Graph @a G.
     * @param G The graph.
     * @param tol The tolerance used for the approximation
     */
    CommuteTimeDistance(const Graph& G, double tol = 0.1);

    /**
     * Destructor.
     */
    ~CommuteTimeDistance() override = default;

    /**
     * Computes ECTD exactly.
     */
    void run() override;

    /**
     * Computes approximation by projection.
     */
    void runApproximation();

    /**
     * Computes approximation by projection, in parallel.
     */
    void runParallelApproximation();

    /**
     * @return The elapsed time to setup the solver in milliseconds.
     */
    uint64_t getSetupTime() const;


    /**
     * Returns the commute time distance between node @a u and node @a v.
     * This method does not need the initial preprocessing (no need to call the run() method).
     * @return commute time distance between the two nodes.
     */

    double runSinglePair(node u, node v);

    /**
     * Returns the commute time distance between node @a u and node @a v.
     * @return commute time distance between the two nodes. Needs to call run() or runApproximation() first.
     */
    double distance(node u, node v);

    double runSingleSource(node u);

protected:
    const Graph* G;
    double tol;
    Lamg<CSRMatrix> lamg;
    uint64_t setupTime;
    std::vector<std::vector<double>> distances;
    std::vector<Vector> solutions;
    bool exactly;
    count k;
};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_COMMUTE_TIME_DISTANCE_HPP_
