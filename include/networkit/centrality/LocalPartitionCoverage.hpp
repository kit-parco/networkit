// no-networkit-format
#ifndef NETWORKIT_CENTRALITY_LOCAL_PARTITION_COVERAGE_HPP_
#define NETWORKIT_CENTRALITY_LOCAL_PARTITION_COVERAGE_HPP_


#include <networkit/centrality/Centrality.hpp>
#include <networkit/structures/Partition.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * The local partition coverage is the amount of neighbors of a node u that are in the same partition as u.
 */
class LocalPartitionCoverage : public Centrality {
public:
    /**
     * Construct the local partition coverage instance. The running time of the run() method is O(m), where m is the number of edges in the graph.
     *
     * @param G The graph to use
     * @param P The partition for which the coverage shall be calculated.
     */
    LocalPartitionCoverage(const Graph& G, const Partition &P);

    /**
     * Computes local partition coverage on the graph passed in constructor.
     * This method runs in parallel.
     */
    void run() override;

    /**
     * Get the maximum value (1.0)
     *
     * @return 1.0
     */
    double maximum() override;

    /**
     * This algorithm is parallel.
     * @return true
     */
    bool TLX_DEPRECATED(isParallel() const override);

    /**
     * The name of this algorithm.
     * @return "Local partition coverage"
     */
    std::string TLX_DEPRECATED(toString() const override);

protected:
    const Partition& P;
};

}

#endif // NETWORKIT_CENTRALITY_LOCAL_PARTITION_COVERAGE_HPP_
