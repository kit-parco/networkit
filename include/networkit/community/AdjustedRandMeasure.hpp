// no-networkit-format
#ifndef NETWORKIT_COMMUNITY_ADJUSTED_RAND_MEASURE_HPP_
#define NETWORKIT_COMMUNITY_ADJUSTED_RAND_MEASURE_HPP_

#include <networkit/community/DissimilarityMeasure.hpp>

namespace NetworKit {

/**
 * The adjusted rand dissimilarity measure as proposed by Huber and Arabie in "Comparing partitions" (http://link.springer.com/article/10.1007/BF01908075)
 */
class AdjustedRandMeasure final : public DissimilarityMeasure {
public:
    /**
     * Get the adjust rand dissimilarity. Runs in O(n log(n)).
     *
     * Note that the dissimilarity can be larger than 1 if the partitions are more different than expected in the random model.
     *
     * @param G    The graph on which the partitions shall be compared
     * @param zeta The first partiton
     * @param eta  The second partition
     * @return The adjusted rand dissimilarity.
     */
    double getDissimilarity(const Graph &G, const NetworKit::Partition &zeta, const NetworKit::Partition &eta) override;
};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_ADJUSTED_RAND_MEASURE_HPP_
