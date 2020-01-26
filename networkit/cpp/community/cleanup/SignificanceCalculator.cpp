/*
 * SignificanceCalculator.cpp
 *
 * Created: 2019-09-13
 * Author: Armin Wiebigke
 */

#include <algorithm>
#include <tuple>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/community/cleanup/SignificanceCalculator.hpp>

namespace NetworKit {

SignificanceCalculator::SignificanceCalculator(const StochasticDistributionCalculator &stochasticDistribution)
        : dist(stochasticDistribution), rng(Aux::Random::integer()),
          random_distribution(-0.5, 0.5) {
}

double
SignificanceCalculator::rScore(count k, count kIn, count cOut, count extStubs) {
    assert(kIn <= cOut);
    count openStubs = extStubs + cOut;
    assert(openStubs >= k);
    ensureMaxValue(openStubs);

    double exactProbability = 0, rightCumulative = 0;
    bool lowSelfLoopProbability = k * k < openStubs;
    if (lowSelfLoopProbability) {
        // Use approximation with hypergeometric distribution
        exactProbability = dist.hypergeometricDistribution(openStubs, cOut, k, kIn);
        rightCumulative = dist.rightCumulativeHypergeometric(openStubs, cOut, k, kIn + 1);
    } else {
        // Calculate the probability using the original distribution
        std::tie(exactProbability, rightCumulative) =
                dist.rightCumulativeStochastic(k, kIn, cOut, extStubs);
        rightCumulative -= exactProbability;
    }

    double bootRandomness = random_distribution(rng) * 1e-6;
    double bootInterval = (0.5 + bootRandomness) *
            dist.hypergeometricDistribution(openStubs, cOut, k, kIn);
//	double bootInterval = (0.5 + bootRandomness) * exactProb; // TODO: Use this instead of hypergeom.
    double score = rightCumulative + bootInterval;
//	double score = rightCum + exactProb; // TODO: Use bootInterval or not?
//	assert(score <= 1.001);

    score = std::min(score, 1.);
    score = std::max(score, 1e-100);
    return score;
}

double SignificanceCalculator::orderStatistic(double rScore, count externalNodes, count pos) {
    assert(pos > 0);
    ensureMaxValue(externalNodes);
    return dist.rightCumulativeBinomial(rScore, externalNodes, pos);
}

void SignificanceCalculator::ensureMaxValue(count maxValue) {
    if (dist.maxValue() < maxValue)
        throw std::runtime_error("Maximum value of the distribution is not high enough.");
}

} /* namespace NetworKit */
