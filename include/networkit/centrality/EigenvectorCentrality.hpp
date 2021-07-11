// no-networkit-format
/*
 * EigenvectorCentrality.hpp
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

#ifndef NETWORKIT_CENTRALITY_EIGENVECTOR_CENTRALITY_HPP_
#define NETWORKIT_CENTRALITY_EIGENVECTOR_CENTRALITY_HPP_

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 * Computes the leading eigenvector of the graph's adjacency matrix (normalized in 2-norm).
 * Interpreted as eigenvector centrality score.
 */
class EigenvectorCentrality final: public Centrality {
    const double tol; // error tolerance

public:
    /**
     * Constructs an EigenvectorCentrality object for the given Graph @a G. @a tol defines the tolerance for convergence.
     *
     * @param[in] G The graph.
     * @param[in] tol The tolerance for convergence.
     * TODO running time
     */
    EigenvectorCentrality(const Graph& G, double tol = 1e-8);

    /**
     * Computes eigenvector centrality on the graph passed in constructor.
     */
    void run() override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_CENTRALITY_EIGENVECTOR_CENTRALITY_HPP_
