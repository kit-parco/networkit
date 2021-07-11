// no-networkit-format
#ifndef NETWORKIT_COMMUNITY_PARTITION_INTERSECTION_HPP_
#define NETWORKIT_COMMUNITY_PARTITION_INTERSECTION_HPP_

#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * Class for calculating the intersection of two partitions, i.e. the clustering with the fewest clusters
 * such that each cluster is a subset of a cluster in both partitions.
 */
class PartitionIntersection final {
    public:
        /**
         * Calculate the intersection of two partitions @a zeta and @a eta
         * @param zeta The first partition
         * @param eta The second partition
         * @return The intersection of @a zeta and @a eta
         */
        Partition calculate(const Partition &zeta, const Partition &eta);
};

} /* namespace NetworKit */

#endif // NETWORKIT_COMMUNITY_PARTITION_INTERSECTION_HPP_
