// no-networkit-format
/*
 * LaplacianCentrality.h
 *
 *  Created on: 08.03.2018
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_CENTRALITY_LAPLACIAN_CENTRALITY_HPP_
#define NETWORKIT_CENTRALITY_LAPLACIAN_CENTRALITY_HPP_

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 * Computes the Laplacian centrality of the graph.
 *
 * The implementation is a simplification of the original algorithm proposed by Qi et al. in
 * "Laplacian centrality: A new centrality measure for weighted networks".
 *
 * See https://dl.acm.org/citation.cfm?id=2181343.2181780 for details.
 */
class LaplacianCentrality: public Centrality {
public:
    /**
     * Constructs a LaplacianCentrality object for the given Graph @a G.
     *
     * @param G The graph.
     * @param normalized Whether scores should be normalized by the energy of the full graph.
     */
    LaplacianCentrality(const Graph& G, bool normalized = false);

    /**
     * Computes the Laplacian centrality on the graph passed in the constructor.
     *
     * See https://dl.acm.org/citation.cfm?id=2181343.2181780 for more details about
     * Laplacian centrality.
     */
    void run() override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_CENTRALITY_LAPLACIAN_CENTRALITY_HPP_
