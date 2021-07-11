#ifndef NETWORKIT_SCD_CLIQUE_DETECT_HPP_
#define NETWORKIT_SCD_CLIQUE_DETECT_HPP_

#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {

/**
 * The CliqueDetect algorithm. It finds the largest clique in the seed node's neighborhood.
 *
 * The algorithm can handle weighted graphs. There, the clique with the highest sum of internal edge
 * weights is returned. This sum includes edge weights to the seed node(s) to ensure that cliques
 * that are well-connected to the seed node(s) are preferred.
 *
 * See also: Hamann, M.; Röhrs, E.; Wagner, D. Local Community Detection Based on Small Cliques.
 * Algorithms 2017, 10, 90. https://doi.org/10.3390/a10030090
 */
class CliqueDetect : public SelectiveCommunityDetector {

public:
    /**
     * Construct a Cliquedetect object.
     *
     * @param[in] G The graph to detect communities on
     */
    CliqueDetect(const Graph &g);

    /**
     * Expands a single seed node/vertex into a maximal clique.
     *
     * @param[in] s the seed node
     * @return A community of the seed node
     */
    std::set<node> expandOneCommunity(node seed) override;

    /**
     * Detect a single clique for the given seed nodes.
     *
     * The resulting community is a clique iff the seeds form a clique.
     * Otherwise, only the added nodes form a clique that is fully connected
     * to the seed nodes.
     *
     * @param seeds The seeds for the community.
     * @return The found community as set of nodes.
     */
    std::set<node> expandOneCommunity(const std::set<node> &seeds) override;

protected:
    std::vector<node> getMaximumWeightClique(const std::vector<node> &nodes,
                                             const std::vector<edgeweight> &seedToNodeWeight) const;
};

} // namespace NetworKit

#endif // NETWORKIT_SCD_CLIQUE_DETECT_HPP_
