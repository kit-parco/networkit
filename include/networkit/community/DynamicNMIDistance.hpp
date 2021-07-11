// no-networkit-format
/*
 * DynamicNMIDistance.hpp
 *
 *  Created on: Jun 26, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_COMMUNITY_DYNAMIC_NMI_DISTANCE_HPP_
#define NETWORKIT_COMMUNITY_DYNAMIC_NMI_DISTANCE_HPP_

#include <networkit/community/DissimilarityMeasure.hpp>
#include <networkit/community/NMIDistance.hpp>

namespace NetworKit {

using Matrix = std::vector<std::vector<count>>;

/**
 * @ingroup community
 */
class DynamicNMIDistance final: public DissimilarityMeasure {
public:

    /**
     * Computes NMI between two clusterings that belong to two different graphs.
     * @a newGraph has evolved from oldGraph, which is only given implicitly via
     * @a oldClustering. NMI is only applied to nodes that belong to the intersection
     * of oldGraph and @a newGraph. Nodes of oldGraph not existing in @newGraph are
     * marked by the entry none in @a newClustering.
     */
    double getDissimilarity(const Graph& newGraph, const Partition& oldClustering, const Partition& newClustering) override;

    void combineValues(double H_sum, double MI, double& NMI, double& NMID) const;
    void sanityCheck(double& NMI, double& NMID) const;

    double entropy(const Partition& clustering, count n, std::vector<double> probs);

    bool isInBoth(node u, const Partition& oldClustering, const Partition& newClustering);

    Matrix confusionMatrix(const Graph& G, const Partition& zeta, const Partition& eta);
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_DYNAMIC_NMI_DISTANCE_HPP_
