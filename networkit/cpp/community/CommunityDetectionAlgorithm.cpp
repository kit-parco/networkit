// no-networkit-format
/*
 * Clusterer.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt
 */

#include <networkit/community/CommunityDetectionAlgorithm.hpp>

namespace NetworKit {

CommunityDetectionAlgorithm::CommunityDetectionAlgorithm(const Graph& G) : Algorithm(), G(&G), result(0) {
    // currently our community detection methods are not defined on directed graphs
    if (G.isDirected()) {
        throw std::runtime_error("This community detection method is undefined on directed graphs");
    }
}

CommunityDetectionAlgorithm::CommunityDetectionAlgorithm(const Graph &G, Partition baseClustering)
    : Algorithm(), G(&G), result(std::move(baseClustering)) {}

Partition CommunityDetectionAlgorithm::getPartition() {
    if(!hasRun) {
        throw std::runtime_error("Call run()-function first.");
    }
    return result;
}

std::string CommunityDetectionAlgorithm::toString() const {
    return "TODO: string representation of clusterer";
}

} /* namespace NetworKit */
