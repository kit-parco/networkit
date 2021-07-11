/*
 * LocalMaxMatcher.hpp
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_MATCHING_LOCAL_MAX_MATCHER_HPP_
#define NETWORKIT_MATCHING_LOCAL_MAX_MATCHER_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/matching/Matcher.hpp>

namespace NetworKit {

/**
 * @ingroup matching
 * LocalMax matching similar to the one described in the EuroPar13 paper
 * by the Sanders group (Birn, Osipov, Sanders, Schulz, Sitchinava)
 */
class LocalMaxMatcher final : public Matcher {
public:
    /**
     * @param G Undirected graph with for which the matching is computed.
     */
    LocalMaxMatcher(const Graph &G);

    /**
     * Runs the algorithm.
     */
    void run() override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_MATCHING_LOCAL_MAX_MATCHER_HPP_
