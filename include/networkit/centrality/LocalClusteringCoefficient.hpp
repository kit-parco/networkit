// no-networkit-format
/*
 * LocalClusteringCoefficient.h
 *
 *  Created on: 31.03.2015
 *      Author: maxv
 */

#ifndef NETWORKIT_CENTRALITY_LOCAL_CLUSTERING_COEFFICIENT_HPP_
#define NETWORKIT_CENTRALITY_LOCAL_CLUSTERING_COEFFICIENT_HPP_

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class LocalClusteringCoefficient: public Centrality {
public:
    /**
     * Constructs the LocalClusteringCoefficient class for the given Graph @a G. If the local clustering coefficient scores should be normalized,
     * then set @a normalized to <code>true</code>. The graph may not contain self-loops.
     *
     * There are two algorithms available. The trivial (parallel) algorithm needs only a small amount of additional memory.
     * The turbo mode adds a (sequential, but fast) pre-processing step using ideas from [0]. This reduces the running time
     * significantly for most graphs. However, the turbo mode needs O(m) additional memory. In practice this should be a bit
     * less than half of the memory that is needed for the graph itself. The turbo mode is particularly effective for graphs
     * with nodes of very high degree and a very skewed degree distribution.
     *
     * [0] Triangle Listing Algorithms: Back from the Diversion
     * Mark Ortmann and Ulrik Brandes                                                                          *
     * 2014 Proceedings of the Sixteenth Workshop on Algorithm Engineering and Experiments (ALENEX). 2014, 1-8
     *
     * @param G The graph.
     * @param turbo If the turbo mode shall be activated.
     * TODO running time
     */
    LocalClusteringCoefficient(const Graph &G, bool turbo = false);



    /**
     * Computes the local clustering coefficient on the graph passed in constructor.
     */
    void run() override;


    /**
    * Get the theoretical maximum of centrality score in the given graph.
    *
    * @return The maximum centrality score.
    */
    double maximum() override;

protected:
    bool turbo;

};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_LOCAL_CLUSTERING_COEFFICIENT_HPP_
