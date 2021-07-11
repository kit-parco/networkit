// no-networkit-format
/*
 * ApproxBetweenness.hpp
 *
 *  Created on: 09.04.2014
 *      Author: cls
 */

#ifndef NETWORKIT_CENTRALITY_APPROX_BETWEENNESS_HPP_
#define NETWORKIT_CENTRALITY_APPROX_BETWEENNESS_HPP_

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {


/**
 * @ingroup centrality
 * Approximation of betweenness centrality according to algorithm described in
 * Matteo Riondato and Evgenios M. Kornaropoulos: Fast Approximation of Betweenness Centrality through Sampling
 */
class ApproxBetweenness final: public Centrality {

public:

    /**
     * The algorithm approximates the betweenness of all vertices so that the scores are
     * within an additive error @a epsilon with probability at least (1- @a delta).
     * The values are normalized by default. The run() method takes O(m) time per sample, where  m is
     * the number of edges of the graph. The number of samples is proportional to universalConstant/epsilon^2.
     * Although this algorithm has a theoretical guarantee, the algorithm implemented in Estimate Betweenness usually performs better in practice
     * Therefore, we recommend to use EstimateBetweenness if no theoretical guarantee is needed.
     *
     * @param	G			the graph
     * @param	epsilon		maximum additive error
     * @param	delta		probability that the values are within the error guarantee
     * @param   universalConstant   the universal constant to be used in
     * computing the sample size. It is 1 by default. Some references suggest
     * using 0.5, but there is no guarantee in this case.
     */
    ApproxBetweenness(const Graph &G, double epsilon = 0.01, double delta = 0.1,
                      double universalConstant = 1.0);

    /**
     * Computes betweenness approximation on the graph passed in constructor.
     */
    void run() override;

    /**
     * @return number of samples taken in last run
     */
    count numberOfSamples() const;


private:

    const double epsilon;
    const double delta;
    count r; // number of samples taken in last run
    double universalConstant;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_APPROX_BETWEENNESS_HPP_
