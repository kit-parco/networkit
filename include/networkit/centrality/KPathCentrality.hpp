// no-networkit-format
/*
 * KPathCentrality.h
 *
 *  Created on: 05.10.2014
 *      Author: nemes
 */

#ifndef NETWORKIT_CENTRALITY_K_PATH_CENTRALITY_HPP_
#define NETWORKIT_CENTRALITY_K_PATH_CENTRALITY_HPP_

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class KPathCentrality: public Centrality {
public:

    /*
     * the maximum length of paths
     * default value ln(n+m)
     */
    count k;
    /*
     * value in interval [-0.5, 0.5]
     * tradeoff between runtime and precision
     * -0.5: maximum precision, maximum runtime
     *  0.5: lowest precision, lowest runtime
     * default value 0.2
     */
    double alpha;

    /**
     * Constructs the K-Path Centrality class for the given Graph @a G.
     *
     * @param G The graph.
     * @param alpha tradeoff between precision and runtime.
     * @param k maximum length of paths.
     * TODO running times
     */
    KPathCentrality(const Graph& G, double alpha=0.2, count k=0);

    /**
     * Computes k-path centrality on the graph passed in constructor.
     */
    void run() override;

};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_K_PATH_CENTRALITY_HPP_
