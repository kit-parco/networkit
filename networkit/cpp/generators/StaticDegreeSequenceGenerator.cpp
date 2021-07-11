/*
 * StaticDegreeSequenceGenerator.cpp
 *
 *  Created on: 24.02.2014
 *      Author: Henning
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/generators/StaticDegreeSequenceGenerator.hpp>

#include <tlx/define.hpp>

namespace NetworKit {

StaticDegreeSequenceGenerator::StaticDegreeSequenceGenerator(const std::vector<count> &sequence)
    : seq(sequence), realizable(UNKNOWN) {}

bool StaticDegreeSequenceGenerator::getRealizable() const {
    return realizable;
}

bool StaticDegreeSequenceGenerator::isRealizable() {
    DEBUG("check if sequence is realizable");
    count n = seq.size();

    // First inequality
    count deg_sum = 0;
    for (count i = 0; i < n; ++i) {
        if (TLX_UNLIKELY(seq[i] >= n)) {
            realizable = NO;
            DEBUG("not realizable: ", seq[i], ", n: ", n);
            return false;
        }
        deg_sum += seq[i];
    }

    if (deg_sum % 2 != 0) {
        DEBUG("not realizable, degree sum not even!");
        realizable = NO;
        return false;
    }

    /**
     * Second inequality
     * We now check that for all 0 <= j < n the following inequality holds
     *   sum(d[i] for 0 <= i <= j) <= sum( min(j+1, d[i]) for j < i < n),
     * where d is the sorted degree sequence.
     *
     * To avoid the quadratic runtime of the naive implementation above, we search for each
     * j the smallest index k with d[k] < j+1. Then the RHS becomes:
     *     (j+1)*min(0, k-j-1)       // all cases where the min-term defaulted to (j+1)
     *   + sum(d[i] for k <= i < n)  // "true" suffix sum.
     *
     * Observe that the suffix sum only depends on k and not j; so we can precompute it once
     * to avoid recomputation.
     *
     * Together, both tricks reduce the runtime complexity from Theta(n^2) to O(n log n).
     */
    std::vector<count> partialSeqSum(n + 1);
    std::copy(seq.cbegin(), seq.cend(), partialSeqSum.begin());
    Aux::Parallel::sort(partialSeqSum.begin(), partialSeqSum.end(), std::greater<count>());
    for (size_t i = n - 1;
         i--;) { // not using std::partial_sum as unclear whether input/output may be identical
        partialSeqSum[i] += partialSeqSum[i + 1];
    }

    auto degreeOf = [&](size_t i) {
        assert(i < n);
        return partialSeqSum[i] - partialSeqSum[i + 1];
    };

    deg_sum = 0;
    for (count j = 0; j < n; ++j) {
        deg_sum += degreeOf(j);
        count min_deg_sum = 0;

        size_t sumFrom = j + 1;
        if (sumFrom < n && degreeOf(sumFrom) >= j + 1) {
            // find the first element right of j that has a value less or equal to j
            const auto it =
                std::lower_bound(partialSeqSum.data() + sumFrom, partialSeqSum.data() + n, j,
                                 [](const count &x, const count j) { return x - *(&x + 1) > j; });
            sumFrom = std::distance(partialSeqSum.data(), it);
            min_deg_sum += (j + 1) * (sumFrom - j - 1);
        }

        if (sumFrom != n)
            min_deg_sum += partialSeqSum[sumFrom];

        if (TLX_UNLIKELY(deg_sum > (j + 1) * j + min_deg_sum)) {
            DEBUG("not realizable");
            realizable = NO;
            return false;
        }
    }

    realizable = YES;
    return true;
}

} /* namespace NetworKit */
