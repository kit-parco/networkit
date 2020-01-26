/*
 * PLP.hpp
 *
 *  Created on: 07.12.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_COMMUNITY_PLP_HPP_
#define NETWORKIT_COMMUNITY_PLP_HPP_

#include <networkit/community/CommunityDetectionAlgorithm.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/community/ClusteringFunctionFactory.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * As described in Ovelgoenne et al: An Ensemble Learning Strategy for Graph Clustering
 * Raghavan et al. proposed a label propagation algorithm for graph clustering.
 * This algorithm initializes every vertex of a graph with a unique label. Then, in iterative
 * sweeps over the set of vertices the vertex labels are updated. A vertex gets the label
 * that the maximum number of its neighbors have. The procedure is stopped when every vertex
 * has the label that at least half of its neighbors have.
 *
 */
class PLP final : public CommunityDetectionAlgorithm {

private:

    count updateThreshold = 0;
    count maxIterations;
    count nIterations = 0; //!< number of iterations in last run
    std::vector<count> timing; //!< running times for each iteration

public:

    /**
     * Constructor to the label propagation community detection algorithm.
     *
     * @param[in] G input graph
     * @param[in] theta optional; updateThreshold: number of nodes that have to be changed in each iteration so that a new iteration starts.
     * @param[in] maxIterations optional; the maximum number of iterations
     */
    PLP(const Graph& G, count theta = none, count maxIterations=none);

    /**
     * Constructor to the label propagation community detection algorithm, starting from a given
     * clustering.
     *
     * @param[in] G input graph
     * @param[in] baseClustering optional; the algorithm will start from the given clustering.
     * @param[in] theta optional; updateThreshold: number of nodes that have to be changed in each iteration so that a new iteration starts.
     * @param[in] maxIterations optional; the maximum number of iterations
     */
    PLP(const Graph& G, const Partition& baseClustering, count theta = none,
            count maxIterations=none);

    /**
     * Run the label propagation clustering algorithm.
     */
    void run() override;

    /**
     * @return String representation of algorithm and parameters.
     */
    std::string toString() const override;

    /**
     * The algorithm runs until a number of nodes less than
     * the threshold is updated.
     *
     * @param th The threshold.
    */
    void setUpdateThreshold(count th);

    /**
    * Get number of iterations in last run.
    *
    * @return The number of iterations.
    */
    count numberOfIterations();

    /**
    * Get list of running times for each iteration.
    *
    * @return The list of running times in milliseconds
    */
    std::vector<count> getTiming();

};

class PLPFactory : public ClusteringFunctionFactory {
public:
    explicit PLPFactory(count theta = none, count maxIterations=none);

    ClusteringFunction getFunction() const override;

private:
    count theta;
    count maxIterations;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_PLP_HPP_
