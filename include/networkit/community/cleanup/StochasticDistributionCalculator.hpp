/*
 * StochasticDistributionCalculator.hpp
 *
 * Created: 2019-09-13
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_COMMUNITY_CLEANUP_STOCHASTIC_DISTRIBUTION_CALCULATOR_HPP_
#define NETWORKIT_COMMUNITY_CLEANUP_STOCHASTIC_DISTRIBUTION_CALCULATOR_HPP_

#include <vector>
#include <cassert>

#include <networkit/Globals.hpp>

namespace NetworKit {

/**
 * Efficient calculation of stochastic distributions.
 * This class stores the logarithm sum of all numbers up to a given maximum. This allows the
 * calculation of binomial coefficients in constant time and in extension a fast calculation
 * of stochastic distributions. However, the input values of the distributions must not exceed
 * the set maximum value.
 */
class StochasticDistributionCalculator {
public:
    /**
     * @param maxValue maximum value that can be used as an input value of the distributions
     */
    explicit StochasticDistributionCalculator(index maxValue);

    // Deleted to prevent accidentally copy
    StochasticDistributionCalculator(StochasticDistributionCalculator const &) = delete;

    /**
     * Set the maximum possible value that can be used as an input for the distributions.
     * @param maxValue new maximum value
     */
    void setMaxValue(count maxValue);

    /**
     * Increase the maximum possible value that can be used as an input for the distributions.
     * If the maximum value is already large enough, this does nothing.
     * @param maxValue new maximum value
     */
    void increaseMaxValueTo(count maxValue);

    /**
     * Returns the maximum value that can be used as an input value.
     */
    count maxValue() const;

    /**
     * Calculate the binomial coefficient "n choose k". Returns a floating point number that may
     * slightly differ from the exact integer result.
     */
    double binomialCoefficient(count n, count k) const;

    /**
     * Calculate the binomial distribution for k success.
     * @param p probability of a success
     * @param n number of trials
     * @param k number of successful trials
     * @return probability of k successes
     */
    double binomialDistribution(double p, count n, count k) const;

    /**
     * Calculate the cumulative binomial distribution that there are k or more success.
     * @param p probability of a success
     * @param n number of trials
     * @param k number of successful trials
     * @return probability of k or more successes
     */
    double rightCumulativeBinomial(double p, count n, count k) const;

    /**
     * Calculate the cumulative binomial distribution that there are k or less success.
     * @param p probability of a success
     * @param n number of trials
     * @param k number of successful trials
     * @return probability of k or less successes
     */
    double leftCumulativeBinomial(double p, count n, count k) const;

    /**
     * Calculate a hypergeometric distribution
     * @param N number of elements
     * @param K number of successful elements
     * @param n number of trials
     * @param k number of successes
     * @return probability of k successes
     */
    double hypergeometricDistribution(count N, count K, count n, count k) const;

    /**
     * Calculate the cumulative hypergeometric distribution that there are k or more success.
     * @param N number of elements
     * @param K number of successful elements
     * @param n number of trials
     * @param k number of successes
     * @return probability of k or more successes
     */
    double rightCumulativeHypergeometric(count N, count K, count n, count k) const;

    /**
     * Calculate the cumulative hypergeometric distribution that there are k or less success.
     * @param N number of elements
     * @param K number of successful elements
     * @param n number of trials
     * @param k number of successes
     * @return probability of k or less successes
     */
    double leftCumulativeHypergeometric(count N, count K, count n, count k) const;

    /**
     * Calculate the probability that a node has kIn or more edges to the community, according to
     * the null model, as described in the OSLOM paper.
     * @param kTotal degree of the node
     * @param kIn number of edges between node and community
     * @param cOut number of outgoing stubs from the community
     * @param extStubs number of stubs in the rest of the graph
     * @return a pair (p(x = kIn), p(x >= kIn))
     */
    std::pair<double, double> rightCumulativeStochastic(count kTotal, count kIn, count cOut,
                                                        count extStubs) const;

private:
    std::vector<double> logSum; // logSum[x] = Sum(log(i)), i = 1 to x
    static constexpr double precision = 1e-6;

    /**
     * @return natural logarithm of the binomial coefficient "n choose k"
     */
    inline double logBinomCoeff(count n, count k) const;
};

} /* namespace NetworKit */

#endif // NETWORKIT_COMMUNITY_CLEANUP_STOCHASTIC_DISTRIBUTION_CALCULATOR_HPP_
