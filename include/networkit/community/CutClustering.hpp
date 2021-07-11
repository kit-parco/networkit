// no-networkit-format
#ifndef NETWORKIT_COMMUNITY_CUT_CLUSTERING_HPP_
#define NETWORKIT_COMMUNITY_CUT_CLUSTERING_HPP_

#include <networkit/community/CommunityDetectionAlgorithm.hpp>

namespace NetworKit {

/**
 * Cut clustering algorithm as defined in
 * Flake, Gary William; Tarjan, Robert E.; Tsioutsiouliklis, Kostas. Graph Clustering and Minimum Cut Trees.
 * Internet Mathematics 1 (2003), no. 4, 385--408.
 */
class CutClustering final : public CommunityDetectionAlgorithm {
public:
    /**
     * Initialize cut clustering algorithm with parameter alpha.
     *
     * A value of 0 gives a clustering with one cluster with all nodes,
     * A value that equals to the largest edge weight gives singleton clusters.
     *
     * @param alpha The parameter for the cut clustering
     */
    CutClustering(const Graph& G, edgeweight alpha);

    /**
     * Apply algorithm to graph
     *
     * Warning: due to numerical errors the resulting clusters might not be correct.
     * This implementation uses the Edmonds-Karp algorithm for the cut calculation.
     */
    void run() override;

    /**
     * @return string representation of algorithm and parameters.
     */
    std::string TLX_DEPRECATED(toString() const override);

    /**
     * Get the complete hierarchy with all possible parameter values.
     *
     * Each reported parameter value is the lower bound for the range in which the corresponding clustering is calculated by the cut clustering algorithm.
     *
     * Warning: all reported parameter values are slightly too high in order to avoid wrong clusterings because of numerical inaccuracies.
     * Furthermore the completeness of the hierarchy cannot be guaranteed because of these inaccuracies.
     * This implementation hasn't been optimized for performance.
     *
     * @param G The Graph instance for which the hierarchy shall be calculated
     *
     * @return The hierarchy as map
     */
    static std::map<edgeweight, Partition> getClusterHierarchy(const Graph& G);

private:

    /**
     * Helper function for the recursive clustering hierarchy calculation.
     */
    static void clusterHierarchyRecursion(const Graph &G, edgeweight lower, Partition lowerClusters, edgeweight upper, Partition upperClusters, std::map< edgeweight, Partition > &result);
    edgeweight alpha;
};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_CUT_CLUSTERING_HPP_
