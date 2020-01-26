/*
 * StochasticDistributionCalculator.cpp
 *
 * Created: 2019-09-13
 * Author: Armin Wiebigke
 */

#include <cassert>
#include <cmath>
#include <algorithm>

#include <networkit/community/cleanup/StochasticDistributionCalculator.hpp>

namespace NetworKit {

StochasticDistributionCalculator::StochasticDistributionCalculator(index maxValue) {
    setMaxValue(maxValue);
}

void StochasticDistributionCalculator::setMaxValue(count maxValue) {
    if (maxValue <= logSum.size()) {
        logSum.resize(maxValue + 1);
        return;
    }
    logSum.reserve(maxValue + 1);
    if (logSum.empty())
        logSum.push_back(0.0);
    size_t old_size = logSum.size();
    double currentValue = logSum.back();
    logSum.resize(maxValue + 1);
#pragma omp parallel for schedule(static)
    for (omp_index i = old_size; i <= static_cast<omp_index>(maxValue); ++i) {
        logSum[i] = std::log(i);
    }
    for (index i = old_size; i <= maxValue; ++i) {
        currentValue += logSum[i];
        logSum[i] = currentValue;
    }
    assert(logSum.size() == maxValue + 1);
}

void StochasticDistributionCalculator::increaseMaxValueTo(count maxValue) {
    if (maxValue > logSum.size()) {
        setMaxValue(maxValue);
    }
}

count StochasticDistributionCalculator::maxValue() const {
    return logSum.size() - 1;
}

double StochasticDistributionCalculator::binomialCoefficient(count n, count k) const {
    return std::exp(logBinomCoeff(n, k));
}

double StochasticDistributionCalculator::binomialDistribution(double p, count n, count k) const {
    return std::exp(logBinomCoeff(n, k) + k * std::log(p) + (n - k) * std::log(1 - p));
}

// Calculates the change ratio of the binomial coefficient "n choose k" if we increment k to k+1
// => returns ("n choose k+1") / ("n choose k")
inline double binomialCoeffChangeForIncrement(count n, count k) {
    return (double) (n - k) / (k + 1);
}

// Calculates the change ratio of the binomial coefficient "n choose k" if we decrement k to k-1
// => returns ("n choose k-1") / ("n choose k")
inline double binomialCoeffChangeForDecrement(count n, count k) {
    return (double) k / (n + 1 - k);
}

// Calculates the change ratio of a binomial distribution if k is incremented (k -> k+1) or
// decremented (k -> k-1).
class BinomialChangeRatio {
public:
    BinomialChangeRatio(double p, double n) : n(n) {
        successChangeIncrementing = p / (1 - p);
        successChangeDecrementing = (1 - p) / p;
    }

    double incrementingK(count k) {
        return successChangeIncrementing
               * binomialCoeffChangeForIncrement(n, k);
    }

    double decrementingK(count k) {
        return successChangeDecrementing
               * binomialCoeffChangeForDecrement(n, k);
    }

private:
    count n;
    double successChangeIncrementing;
    double successChangeDecrementing;
};

double StochasticDistributionCalculator::rightCumulativeBinomial(double p, count n, count k) const {
    assert("Error: k > n" && k <= n);
    assert("Error: p < 0" && p >= 0);
    if (k == 0)
        return 1;
    if (p - 1 > -1e-11)
        return 1;
    // If k is smaller than the expected value, then calculating the left cumulative is faster
    if (k < n * p)
        return 1 - leftCumulativeBinomial(p, n, k - 1);

    double startBinom = binomialDistribution(p, n, k);
    if (startBinom <= 1e-40)
        return 0;

    // Sum the probabilities until the sum only changes minimally.
    // To increase the performance, we calculate the change in the probability of the
    // distribution for increasing k instead of calculating each probability directly.
    // At the end, the sum is normalized with the starting probability.
    double currentBinom = 1.;
    double sum = currentBinom;
    BinomialChangeRatio changeRatio(p, n);
    for (count x = k; x < n; ++x) {
        currentBinom *= changeRatio.incrementingK(x);
        sum += currentBinom;
        if (currentBinom < precision * sum)
            break;
    }
    assert(startBinom * sum <= 1.001);
    return startBinom * sum;
}

double StochasticDistributionCalculator::leftCumulativeBinomial(double p, count n, count k) const {
    if (k == n)
        return 1;
    if (p < 1e-11)
        return 1;
    // If k is larger than the expected value, then calculating the right cumulative is faster
    if (k > n * p)
        return 1 - rightCumulativeBinomial(p, n, k + 1);

    double startBinom = binomialDistribution(p, n, k);
    if (startBinom <= 1e-40)
        return 0;

    // Same approach as in rightCumulativeBinomial
    double curBinom = 1.;
    double sum = curBinom;
    BinomialChangeRatio changeRatio(p, n);
    for (count x = k; x > 0; --x) {
        curBinom *= changeRatio.decrementingK(x);
        sum += curBinom;
        if (curBinom < precision * sum)
            break;
    }
    assert(startBinom * sum <= 1.001);
    return startBinom * sum;
}

double
StochasticDistributionCalculator::hypergeometricDistribution(count N, count K, count n, count k) const {
    double logValue = logBinomCoeff(K, k)
                      + logBinomCoeff(N - K, n - k)
                      - logBinomCoeff(N, n);
    return std::exp(logValue);
}


// Calculates the change ratio of a hypergeometric distribution if k is incremented (k -> k+1) or
// decremented (k -> k-1).
class HypergeomChangeRatio {
public:
    HypergeomChangeRatio(count N, count K, count n) : N(N), K(K), n(n) {};

    double incrementingK(count k) {
        return binomialCoeffChangeForIncrement(K, k) *
               binomialCoeffChangeForDecrement(N - K, n - k);
    }

    double decrementingK(count k) {
        return binomialCoeffChangeForDecrement(K, k) *
               binomialCoeffChangeForIncrement(N - K, n - k);
    }

private:
    count N;
    count K;
    count n;
};

double
StochasticDistributionCalculator::rightCumulativeHypergeometric(count N, count K, count n, count k) const {
    assert(K <= N && "Error: K > N");
    assert(this->maxValue() >= N);
    if (k == 0)
        return 1;
    if (k > n || k > K)
        return 0;
    // If k is smaller than the expected value, then calculating the left cumulative is faster
    double expectedValue = (double) K / N * n;
    if (k < expectedValue)
        return (1. - leftCumulativeHypergeometric(N, K, n, k - 1));

    double startProb = hypergeometricDistribution(N, K, n, k);
    if (startProb <= 1e-40)
        return 0;

    // Same approach as in rightCumulativeBinomial
    double curProb = 1.;
    double sum = curProb;
    HypergeomChangeRatio changeRatio(N, K, n);
    for (count x = k; x <= n; ++x) {
        curProb *= changeRatio.incrementingK(x);
        sum += curProb;
        if (curProb < precision * sum)
            break;
    }

    assert(startProb * sum <= 1.001);
    return startProb * sum;
}

double
StochasticDistributionCalculator::leftCumulativeHypergeometric(count N, count K, count n, count k) const {
    assert(N < logSum.size());
    if ((int64_t) N - (int64_t) K - (int64_t) n + (int64_t) k < 0)
        return 0;
    if (k == n)
        return 1;
    // If k is larger than the expected value, then calculating the right cumulative is faster
    double expectedValue = (double) K / N * n;
    if (k > expectedValue)
        return (1. - rightCumulativeHypergeometric(N, K, n, k + 1));

    double startProb = hypergeometricDistribution(N, K, n, k);
    if (startProb <= 1e-40)
        return 0;

    // Same approach as in rightCumulativeBinomial
    double curProb = 1.;
    double sum = curProb;
    HypergeomChangeRatio changeRatio(N, K, n);
    for (count x = k; x > 0; --x) {
        curProb *= changeRatio.decrementingK(x);
        sum += curProb;
        if (curProb < precision * sum)
            break;
    }

    assert(startProb * sum <= 1.001);
    return startProb * sum;
}

// Calculates the change ratio of the statistical significance distribution if k is
// incremented (k -> k+1) or decremented (k -> k-1).
class SignificanceChangeRatio {
public:
    SignificanceChangeRatio(count kTotal, count cOut, count extStubs) : kTotal(kTotal), cOut(cOut) {
        count M = extStubs - kTotal;
        MInEdgesConstPart =
                (M - cOut - kTotal) / 2; // MIn = M - cOut - kTotal + 2*kIn, MInEdges = M / 2
    }

    double incrementingK(count kIn) {
        return 0.5                    // 2^-kIn
               * (kTotal - kIn)       // 1/kOut! (-1 in Factorial)
               / (kIn + 1)            // 1/kIn!  (+1 in Factorial)
               * (cOut - kIn)         // 1/(cOut - kIn)!  (-1 in Factorial)
               / (MInEdgesConstPart + (kIn + 1))     // 1/(MIn/2)!  (+1 in Factorial)
                ;
    }

    double decrementingK(count kIn) {
        return 2.                           // 2^-kIn
               / (kTotal - (kIn - 1))       // 1/kOut! (+1 in Factorial)
               * kIn                        // 1/kIn!  (-1 in Factorial)
               / (cOut - (kIn - 1))         // 1/(cOut - kIn)!  (+1 in Factorial)
               * (MInEdgesConstPart + kIn)  // 1/(MIn/2)!  (-1 in Factorial)
                ;
    }

private:
    count kTotal;
    count cOut;
    count MInEdgesConstPart;
};

// TODO: Better name?
std::pair<double, double>
StochasticDistributionCalculator::rightCumulativeStochastic(count kTotal, count kIn, count cOut,
                                                            count extStubs) const {
    assert(kTotal >= kIn);
    assert(kTotal < logSum.size() && cOut < logSum.size() && extStubs < logSum.size());
    if (kIn > cOut || kIn > kTotal)
        return {0.0, 0.0};

//	// TODO: Use mode for speedup?
//	// Get mode == highest probability
//	// TODO: Is there a better guess for the mode?
//	count mode = (cOut / double(M + kTotal) * kTotal);
//	mode = std::min(mode, cOut); // this mode is underestimated anyway

    // Sum all (non-negligible) probabilities to get the normalization factor.

    // Right cumulative
    double kInProbability = 1.0; // p(kIn) (not normalized)
    double rightCumulativeSum = kInProbability;
    double currentProbability = kInProbability;
    count maxKIn = std::min(kTotal, cOut);
    SignificanceChangeRatio changeRatio(kTotal, cOut, extStubs);
    for (count x = kIn; x < maxKIn; ++x) {
        double ratio = changeRatio.incrementingK(x);
        currentProbability *= ratio;
        rightCumulativeSum += currentProbability;
        if (rightCumulativeSum > 1e256) {
            // currentProbability can become very large if the right cumulative probability is
            // very low. In this case, we normalize all probabilities.
            currentProbability /= 1e256;
            rightCumulativeSum /= 1e256;
            kInProbability /= 1e256;
        }
        assert(ratio != 0.0);
        assert(currentProbability < 1e256);

        double sumChange = currentProbability / rightCumulativeSum;
        if (sumChange < precision)
            break;
    }

    // Left cumulative
    double probabilitySum = rightCumulativeSum;
    currentProbability = kInProbability;
    for (count x = kIn; x > 0; --x) {
        assert(x < 1e100);
        currentProbability *= changeRatio.decrementingK(x);
        probabilitySum += currentProbability;
        if (probabilitySum > 1e256) {
            // currentProbability can become very large if the right cumulative probability is
            // very low. In this case, we normalize all probabilities.
            currentProbability /= 1e256;
            rightCumulativeSum /= 1e256;
            kInProbability /= 1e256;
            probabilitySum /= 1e256;
        }

        double sumChange = currentProbability / probabilitySum;
        if (sumChange < precision)
            break;
    }

    double normalizedKInProbability = kInProbability / probabilitySum;
    double normalizedRightCumulative = rightCumulativeSum / probabilitySum;
    assert(normalizedKInProbability < 1.001);
    assert(normalizedRightCumulative < 1.001);
    return {normalizedKInProbability, normalizedRightCumulative};
}

double StochasticDistributionCalculator::logBinomCoeff(count n, count k) const {
    assert(n >= k);
    return logSum[n] - logSum[n - k] - logSum[k];
}

} /* namespace NetworKit */
