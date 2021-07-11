// no-networkit-format
#ifndef NETWORKIT_COMMUNITY_PARTITION_FRAGMENTATION_HPP_
#define NETWORKIT_COMMUNITY_PARTITION_FRAGMENTATION_HPP_

#include <networkit/community/LocalPartitionEvaluation.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * This measure evaluates how fragmented a partition is. The fragmentation of a single cluster is defined as one minus the
 * number of nodes in its maximum connected components divided by its total number of nodes. Smaller values thus indicate a smaller fragmentation.
 */
class PartitionFragmentation final : public LocalPartitionEvaluation {
public:
    using LocalPartitionEvaluation::LocalPartitionEvaluation;

    /**
     * Execute the algorithm. The algorithm is not parallel.
     */
    void run() override;

    /**
     * @return true - smaller values are better than larger values.
     */
    bool isSmallBetter() const override { return true; };

    /**
     * @return false - only minor parts of this implementation are parallel.
     */
    bool TLX_DEPRECATED(isParallel() const override) { return false; };
};

}

#endif // NETWORKIT_COMMUNITY_PARTITION_FRAGMENTATION_HPP_
