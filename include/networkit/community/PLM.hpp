/*
 * PLM.hpp
 *
 *  Created on: 20.11.2013
 *      Author: cls
 */

#ifndef NETWORKIT_COMMUNITY_PLM_HPP_
#define NETWORKIT_COMMUNITY_PLM_HPP_

#include <networkit/community/CommunityDetectionAlgorithm.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Parallel Louvain Method - a multi-level modularity maximizer.
 */
class PLM final : public CommunityDetectionAlgorithm {

public:
    /**
     * @param[in] G input graph
     * @param[in] zeta initial partition
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
     *
     */
    PLM(const Graph& G, bool refine=false, double gamma = 1.0, std::string par="balanced", count maxIter=32, bool turbo = true, bool recurse = true, Partition zeta=Partition(0));

    PLM(const Graph& G, const PLM& other);

    /**
     * Get string representation.
     *
     * @return String representation of this algorithm.
     */
    std::string TLX_DEPRECATED(toString() const override);

    /**
     * Detect communities.
     */
    void run() override;

    static std::pair<Graph, std::vector<node>> coarsen(const Graph& G, const Partition& zeta);

    static Partition prolong(const Graph& Gcoarse, const Partition& zetaCoarse, const Graph& Gfine, std::vector<node> nodeToMetaNode);

    /**
     * Returns fine-grained running time measurements for algorithm engineering purposes.
     */
    std::map<std::string, std::vector<count> > getTiming();

private:

    Partition zeta;
    std::string parallelism;
    bool refine;
    double gamma = 1.0;
    count maxIter;
    bool turbo;
    bool recurse;
    std::map<std::string, std::vector<count>> timing;  // fine-grained running time measurement
};

} /* namespace NetworKit */

#endif // NETWORKIT_COMMUNITY_PLM_HPP_
