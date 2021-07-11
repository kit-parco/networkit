// no-networkit-format
/*
 * ParallelAgglomerativeClusterer.hpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt, Henning Meyerhenke
 */

#ifndef NETWORKIT_COMMUNITY_PARALLEL_AGGLOMERATIVE_CLUSTERER_HPP_
#define NETWORKIT_COMMUNITY_PARALLEL_AGGLOMERATIVE_CLUSTERER_HPP_

#include <networkit/community/CommunityDetectionAlgorithm.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * A parallel agglomerative community detection algorithm, maximizing modularity.
 */
class ParallelAgglomerativeClusterer final: public CommunityDetectionAlgorithm {

public:
    /**
     * Constructor to the parallel agglomerative clusterer.
     *
     * @param[in] G input graph
     */
    ParallelAgglomerativeClusterer(const Graph& G);

    /**
     * Detect communities.
     */
    void run() override;

    std::string TLX_DEPRECATED(toString() const override);
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_PARALLEL_AGGLOMERATIVE_CLUSTERER_HPP_
