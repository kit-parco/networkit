// no-networkit-format
/*
 * SpanningEdgeCentrality.h
 *
 *  Created on: 29.07.2015
 *      Author: henningm
 */

#ifndef NETWORKIT_CENTRALITY_SPANNING_EDGE_CENTRALITY_HPP_
#define NETWORKIT_CENTRALITY_SPANNING_EDGE_CENTRALITY_HPP_

#include <networkit/numerics/LAMG/Lamg.hpp>
#include <networkit/centrality/Centrality.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>


namespace NetworKit {

/**
 * @ingroup centrality
 *
 * SpanningEdgeCentrality edge centrality.
 *
 */
class SpanningEdgeCentrality: public Centrality {
protected:
    double tol;
    Lamg<CSRMatrix> lamg;
    uint64_t setupTime;

public:
    /**
     * Constructs the SpanningEdgeCentrality class for the given Graph @a G.
     * @param G The graph.
     * @param tol constant used for the approximation: with probability at least 1-1/n, the approximated scores are within a factor 1+tol from the exact scores
     */
    SpanningEdgeCentrality(const Graph& G, double tol = 0.1);

    /**
     * Destructor.
     */
    ~SpanningEdgeCentrality() override = default;

    /**
    * Compute spanning edge centrality scores exactly for all edges. This solves a linear system for each edge, so the empirical running time is O(m^2),
    * where m is the number of edges in the graph.
    */
    void run() override;


    /**
     * Compute approximation by JL projection. This solves k linear systems, where k is log(n)/(tol^2). The empirical running time is O(km), where n is the number of nodes
     * and m is the number of edges.
     */
    void runApproximation();

    /**
     * Compute approximation by JL projection in parallel. This solves k linear systems in parallel, where k is log(n)/(tol^2).
     */
    void runParallelApproximation();

    /**
     * Only used by benchmarking. Computes an approximation by projection and solving Laplacian systems.
     * Measures the time needed to compute the approximation and writes the problem vectors to the
     * directory of the graph specified by @a graphPath.
     * @param directory
     * @return Elapsed time in milliseconds.
     */
    uint64_t runApproximationAndWriteVectors(const std::string &graphPath);

    /**
     * @return The elapsed time to setup the solver in milliseconds.
     */
    uint64_t getSetupTime() const;
    /**
     * Compute value for one edge only. This requires a single linear system, so the empricial running time is O(m).
     * @param[in] u Endpoint of edge.
     * @param[in] v Endpoint of edge.
     */
    double runForEdge(node u, node v);

};

} /* namespace NetworKit */


#endif // NETWORKIT_CENTRALITY_SPANNING_EDGE_CENTRALITY_HPP_
