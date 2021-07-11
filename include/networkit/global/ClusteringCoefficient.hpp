// no-networkit-format
/*
 * ClusteringCoefficient.hpp
 *
 *  Created on: 08.04.2013
 *      Author: Lukas Barth, David Weiss
 */

#ifndef NETWORKIT_GLOBAL_CLUSTERING_COEFFICIENT_HPP_
#define NETWORKIT_GLOBAL_CLUSTERING_COEFFICIENT_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup global
 */
class ClusteringCoefficient final {

public:
    /**
     * DEPRECATED: use centrality.LocalClusteringCoefficient and take average
     *
     * This calculates the average local clustering coefficient of graph @a G.
     *
     * @param G The graph (may not contain self-loops).
     * @note $$c(G) := \\frac{1}{n} \\sum_{u \\in V} c(u)$$
     * where $c(u) := \\frac{2 \\cdot |E(N(u))| }{\\deg(u) \\cdot ( \\deg(u) - 1)}$
     */
    static double avgLocal(Graph& G, bool turbo = false);
    static double sequentialAvgLocal(const Graph &G);
    static double approxAvgLocal(Graph& G, count trials);

    /**
     * This calculates the global clustering coefficient
     */
    static double exactGlobal(Graph& G);
    static double approxGlobal(Graph& G, count trials);

};

} /* namespace NetworKit */
#endif // NETWORKIT_GLOBAL_CLUSTERING_COEFFICIENT_HPP_
