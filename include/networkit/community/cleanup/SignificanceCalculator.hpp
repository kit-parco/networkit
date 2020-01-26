/*
 * SignificanceCalculator.hpp
 *
 * Created: 2019-09-13
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_COMMUNITY_CLEANUP_SIGNIFICANCE_CALCULATOR_HPP_
#define NETWORKIT_COMMUNITY_CLEANUP_SIGNIFICANCE_CALCULATOR_HPP_

#include <random>
#include <stdexcept>

#include <networkit/Globals.hpp>
#include <networkit/community/cleanup/StochasticDistributionCalculator.hpp>

namespace NetworKit {

/**
 * This class calculates the statistical significance of a node to a community.
 */
class SignificanceCalculator {
public:

    /**
     * Constructor
     * @param dist calculates the stochastic distributions
     */
    explicit SignificanceCalculator(const StochasticDistributionCalculator &dist);

    /**
     * Calculate the r-score
     * @param kIn Number of edges between node and community
     * @param cOut Number of outgoing stubs from the community
     * @param extStubs Number of stubs in the rest of the graph (without the node and the community)
     * @param k Degree of the node
     * @return a pair (r-score, boot interval)
     */
    double rScore(count k, count kIn, count cOut, count extStubs);

    /**
     * Calculate the order statistic (s-score)
     * @param rScore the r-score of the candidate
     * @param externalNodes the number of external nodes
     * @param pos the position of the candidate
     */
    double orderStatistic(double rScore, count externalNodes, count pos);

private:
    const StochasticDistributionCalculator &dist;
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> random_distribution;

    void ensureMaxValue(count maxValue);
};

} /* namespace NetworKit */

#endif // NETWORKIT_COMMUNITY_CLEANUP_SIGNIFICANCE_CALCULATOR_HPP_
