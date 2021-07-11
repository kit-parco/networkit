// no-networkit-format
#ifndef NETWORKIT_COMMUNITY_ISOLATED_INTERPARTITION_CONDUCTANCE_HPP_
#define NETWORKIT_COMMUNITY_ISOLATED_INTERPARTITION_CONDUCTANCE_HPP_

#include <networkit/community/LocalPartitionEvaluation.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * Isolated inter-partition conductance is a measure for how well a partition
 * (communtiy/cluster) is separated from the rest of the graph.
 *
 * The conductance of a partition is defined as the weight of the cut divided
 * by the volume (the sum of the degrees) of the nodes in the partition or the
 * nodes in the rest of the graph, whatever is smaller. Small values thus indicate
 * that the cut is small compared to the volume of the smaller of the separated
 * parts. For the whole partitions usually the maximum or the unweighted average
 * is used.
 *
 * See also Experiments on Density-Constrained Graph Clustering,
 * Robert Gorke, Andrea Kappes and  Dorothea Wagner, JEA 2015:
 * http://dx.doi.org/10.1145/2638551
 */
class IsolatedInterpartitionConductance final : public LocalPartitionEvaluation {
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

    /**
     * Get the name of the algorithm.
     */
    std::string TLX_DEPRECATED(toString() const override) { return "Isolated inter-partition conductance"; };
};

}

#endif // NETWORKIT_COMMUNITY_ISOLATED_INTERPARTITION_CONDUCTANCE_HPP_
