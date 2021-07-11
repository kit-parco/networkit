#ifndef NETWORKIT_COMMUNITY_OVERLAPPING_NMI_DISTANCE_HPP_
#define NETWORKIT_COMMUNITY_OVERLAPPING_NMI_DISTANCE_HPP_

#include <unordered_map>

#include <networkit/auxiliary/HashUtils.hpp>
#include <networkit/community/DissimilarityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * 	Compare two covers using the overlapping normalized mutual information measure. This is a
 * 	dissimilarity measure with a range of [0, 1]. A value of 0 indicates a perfect agreement while
 * 	a 1 indicates complete disagreement.
 *
 * 	For the `MAX` normalization, this is the measure introduced in [NMI13]. Other normalization
 * 	methods result in similar measures.
 *
 * 	Please note that non-overlapping NMIDistance uses another definition of the normalized mutual
 * 	information. See NMIDistance for details on its computation. Both NMIDistance and
 * 	OverlappingNMIDistance can be used with partitions, but produce different values.
 *
 * [NMI13]
 *      McDaid, Aaron F., Derek Greene, and Neil Hurley.
 *      “Normalized Mutual Information to Evaluate Overlapping Community Finding Algorithms.”
 *      ArXiv:1110.2515 [Physics], August 2, 2013. http://arxiv.org/abs/1110.2515.
 */
class OverlappingNMIDistance final : public DissimilarityMeasure {

public:
    enum Normalization { MIN, GEOMETRIC_MEAN, ARITHMETIC_MEAN, MAX, JOINT_ENTROPY };

private:
    Normalization mNormalization{Normalization::MAX};

public:
    OverlappingNMIDistance() = default;

    explicit OverlappingNMIDistance(Normalization normalization) : mNormalization(normalization) {}

    void setNormalization(Normalization normalization) { mNormalization = normalization; }

    double getDissimilarity(const Graph &G, const Partition &zeta, const Partition &eta) override;
    double getDissimilarity(const Graph &G, const Cover &zeta, const Cover &eta) override;

private:
    struct SizesAndIntersections {
        std::vector<count> sizesX;
        std::vector<count> sizesY;
        std::unordered_map<std::pair<index, index>, count, Aux::PairHash> intersectionSizes;
    };

    /**
     * Calculating cluster sizes of covers and intersection sizes between two covers.
     *
     * @param X
     * @param Y
     * @return The cluster and intersection sizes of X and Y.
     */
    static SizesAndIntersections
    calculateClusterAndIntersectionSizes(const Graph &graph, const Cover &X, const Cover &Y);

    /**
     * Calculates partial entropy.
     *
     * $$
     *   h(w, n) = -w \log_2\left(\frac{w}{n}\right)
     * $$
     *
     * @param w  Number of occurrences.
     * @param n  Total number of nodes.
     * @return
     */
    static double h(count w, count n);

    /**
     * Calculates the entropy of a cluster.
     *
     * $$
     *   H(X_i) = h(|\{ u | u \in X_i \}, n) + h(\{ u | u \notin X_i \}, n)
     * $$
     *
     * @param size  Cluster size.
     * @param n  Total number of nodes.
     * @return
     */
    static double entropy(count size, count n);

    /**
     * Calculates the entropy of a clustering.
     *
     * $$
     *   H(X) = \sum_{X_i \in X} H(X_i)
     * $$
     *
     * @param sizesX  The sizes of the clusters in X.
     * @param n  Total number of nodes.
     * @return
     */
    static double entropy(const std::vector<count> &sizesX, count n);

    /**
     * Calculates the adjusted conditional entropy between two clusters.
     *
     * $$
     *   H^*(X_i|Y_j) =
     *   \begin{cases}
     *       H(X_i|Y_j), & \text{if $h(a, n) + h(d, n) \geq h(b, n) + h(c, n)$}\\
     *       H(X_i), & \text{otherwise}
     *   \end{cases}
     * $$
     * where
     * $$
     *   H(X_i|Y_j) &= H(X_i, Y_j) - H(Y_j)\\
     *       &= h(a, n) + h(b, n) + h(c, n) + h(d, n) - h(b + d, n) - h(a + c, n),
     * $$
     * $a = \sum_u \[X_{i,u} = 0 and Y_{j,u} = 0\]$, $b = \sum_u \[X_{i,u} = 0 and Y_{j,u} = 1\]$,
     * $c = \sum_u \[X_{i,u} = 1 and Y_{j,u} = 0\]$ and $d = \sum_u \[X_{i,u} = 1 and Y_{j,u} =
     * 1\]$.
     *
     * @param sizeXi  Size of cluster X_i.
     * @param sizeYj  Size of cluster Y_j.
     * @param intersectionSize  Size of the intersection of X_i and Y_j.
     * @param n  Total number of nodes.
     * @return
     */
    static double adjustedConditionalEntropy(count sizeXi, count sizeYj, count intersectionSize,
                                             count n);

    /**
     * Calculates the entropy of the clustering X conditioned by the clustering Y.
     *
     * $$
     *   H(X|Y) = \sum_{X_i \in X} \min_{Y_j \in Y} H^*(X_i|Y_j)
     * $$
     *
     * @param sizesX
     * @param sizesY
     * @param intersectionSizes
     *      A map which stores the the intersection size between clusters X[i] and Y[j] with indices
     * i and j as keys.
     * @param invertPairIndices
     *      Whether to invert the indices of the intersectionSizes map entries.
     *      This allows reusing the same map for H(X|Y) and H(Y|X).
     * @param n  Total number of nodes.
     * @return
     */
    static double conditionalEntropy(
        const std::vector<count> &sizesX, const std::vector<count> &sizesY,
        const std::unordered_map<std::pair<index, index>, count, Aux::PairHash> &intersectionSizes,
        bool invertPairIndices, count n);

    static void clampBelow(double &value, double lowerBound, const char *format,
                           int printPrecision = 20);

    static void clampAbove(double &value, double upperBound, const char *format,
                           int printPrecision = 20);

    /**
     * Normalize the given mutual information value.
     *
     * @param normalization
     * @param mutualInformation
     * @param entropyX
     * @param entropyY
     * @return
     */
    static double normalize(OverlappingNMIDistance::Normalization normalization,
                            double mutualInformation, double entropyX, double entropyY);
};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_OVERLAPPING_NMI_DISTANCE_HPP_
