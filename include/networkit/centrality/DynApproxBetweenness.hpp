/*
 * DynApproxBetweenness.hpp
 *
 *  Created on: 31.07.2014
 *      Author: ebergamini
 */

// networkit-format
#ifndef NETWORKIT_CENTRALITY_DYN_APPROX_BETWEENNESS_HPP_
#define NETWORKIT_CENTRALITY_DYN_APPROX_BETWEENNESS_HPP_

#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/centrality/Centrality.hpp>
#include <networkit/distance/DynSSSP.hpp>
#include <networkit/dynamics/GraphEvent.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <omp.h>

namespace NetworKit {

// pImpl
class DynApproxBetweennessImpl;
class DynApproxBetweennessImplDir;
/**
 * @ingroup centrality
 * Interface for dynamic approximated betweenness centrality algorithm.
 */
class DynApproxBetweenness : public Centrality, public DynAlgorithm {

public:
    /**
     * The algorithm approximates the betweenness of all vertices so that the scores are
     * within an additive error @a epsilon with probability at least (1- @a delta).
     * The values are normalized by default.
     *
     * @param	G			the graph
     * @param	storePredecessors   keep track of the lists of predecessors?
     * @param	epsilon		maximum additive error
     * @param	delta		probability that the values are within the error guarantee
     * @param	universalConstant	the universal constant to be used in
     * computing the sample size. It is 1 by default. Some references suggest
     * using 0.5, but there is no guarantee in this case.
     */
    DynApproxBetweenness(const Graph &G, double epsilon = 0.01, double delta = 0.1,
                         bool storePredecessors = true, double universalConstant = 1.0);
    /**
     * The algorithm approximates the betweenness of all vertices so that the scores are
     * within an additive error @a epsilon with probability at least (1- @a delta).
     * The values are normalized by default.
     *
     * @param	G			the graph
     * @param	storePredecessors   keep track of the lists of predecessors?
     * @param	epsilon		maximum additive error
     * @param	delta		probability that the values are within the error guarantee
     */
    DynApproxBetweenness(const Graph &G, double epsilon = 0.01, double delta = 0.1,
                         bool storePredecessors = true);

    ~DynApproxBetweenness() override;

    /**
     * Runs the static approximated betweenness centrality algorithm on the initial graph.
     */
    void run() override;

    /**
     * Updates the betweenness centralities after an edge insertions on the graph.
     * Notice: it works only with edge insertions and the graph has to be connected.
     *
     * @param e The edge insertions.
     */
    void update(GraphEvent e) override;

    /**
     * Updates the betweenness centralities after a batch of edge insertions on the graph.
     * Notice: it works only with edge insertions and the graph has to be connected.
     *
     * @param batch The batch of edge insertions.
     */
    void updateBatch(const std::vector<GraphEvent> &batch) override;

    /**
     * Get number of path samples used for last calculation
     */
    count getNumberOfSamples();

private:
    std::unique_ptr<DynApproxBetweennessImpl> impl;
    std::unique_ptr<DynApproxBetweennessImplDir> implDir;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_DYN_APPROX_BETWEENNESS_HPP_
