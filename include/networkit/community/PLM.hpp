/*
 * PLM.hpp
 *
 *  Created on: 20.11.2013
 *      Author: cls
 */

#ifndef NETWORKIT_COMMUNITY_PLM_HPP_
#define NETWORKIT_COMMUNITY_PLM_HPP_

#include <networkit/community/CommunityDetectionAlgorithm.hpp>
#include <networkit/community/ClusteringFunctionFactory.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Parallel Louvain Method - a multi-level modularity maximizer.
 */
class PLM final : public CommunityDetectionAlgorithm {

public:
    /**
     * @param[in] G input graph
     * @param[in] refine add a second move phase to refine the communities
     * @param[in] par parallelization strategy
     * @param[in] gammamulti-resolution modularity parameter:
     *            1.0 -> standard modularity
     *            0.0 -> one community
     *            2m -> singleton communities
     * @param[in] maxIter maximum number of iterations for move phase
     * @param[in] parallelCoarsening use parallel graph coarsening
     * @param[in] turbo faster but uses O(n) additional memory per thread
     * @param[in] recurse use recursive coarsening, see http://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.049902 for some explanations (default: true)
     * @param[in] measure_time Measure the running time of various phases. As timers have some overhead in particular for small graphs, this is disabled by default.
     *
     */
    PLM(const Graph &G, bool refine = false, double gamma = 1.0, std::string par = "balanced", count maxIter = 32,
        bool turbo = true, bool recurse = true, bool measure_time = false);

    PLM(const Graph &G, const PLM &other);

    /**
     * Get string representation.
     *
     * @return String representation of this algorithm.
     */
    std::string toString() const override;

    /**
     * Detect communities.
     */
    void run() override;

    static std::pair<Graph, std::vector<node>> coarsen(const Graph &G, const Partition &zeta, bool parallel);

    static Partition
    prolong(const Graph &Gcoarse, const Partition &zetaCoarse, const Graph &Gfine, std::vector<node> nodeToMetaNode);

    /**
     * Returns fine-grained running time measurements for algorithm engineering purposes.
     */
    std::map<std::string, std::vector<count> > getTiming();

private:

    std::string parallelism;
    bool refine;
    double gamma = 1.0;
    count maxIter;
    bool turbo;
    bool recurse;
    bool measure_time;
    std::map<std::string, std::vector<count> > timing;     // fine-grained running time measurement
};

class PLMFactory : public ClusteringFunctionFactory {
public:
    explicit PLMFactory(bool refine = false, double gamma = 1.0, std::string par = "balanced");

    ClusteringFunction getFunction() const override;

private:
    bool refine;
    double gamma;
    std::string par;
};

} /* namespace NetworKit */

#endif // NETWORKIT_COMMUNITY_PLM_HPP_
