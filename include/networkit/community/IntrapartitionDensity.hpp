// no-networkit-format
#ifndef NETWORKIT_COMMUNITY_INTRAPARTITION_DENSITY_HPP_
#define NETWORKIT_COMMUNITY_INTRAPARTITION_DENSITY_HPP_

#include <networkit/community/LocalPartitionEvaluation.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * The intra-cluster density of a partition is defined as the number of existing edges divided by the number of possible edges.
 * The global value is the sum of all existing intra-cluster edges divided by the sum of all possible intra-cluster edges.
 */
class IntrapartitionDensity final : public LocalPartitionEvaluation {
public:
    using LocalPartitionEvaluation::LocalPartitionEvaluation;

    /**
     * Execute the algorithm. The algorithm is not parallel.
     */
    void run() override;

    /**
     * Get the global intra-cluster density.
     *
     * @return The global intra-cluster density.
     */
    double getGlobal() const { assureFinished(); return globalValue; };

    /**
     * @return false - this algorithm is not parallel.
     */
    bool TLX_DEPRECATED(isParallel() const override) { return false; }

    /**
     * This value should be high in a good clustering.
     * @return false - high values are better than small values.
     */
    bool isSmallBetter() const override { return false; }
private:
    double globalValue;
};

}

#endif // NETWORKIT_COMMUNITY_INTRAPARTITION_DENSITY_HPP_
