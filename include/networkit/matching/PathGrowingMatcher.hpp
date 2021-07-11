// no-networkit-format
/*
 * PathGrowingMatcher.hpp
 *
 *  Created on: Jun 13, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_MATCHING_PATH_GROWING_MATCHER_HPP_
#define NETWORKIT_MATCHING_PATH_GROWING_MATCHER_HPP_

#include <networkit/matching/Matcher.hpp>
#include <networkit/matching/Matching.hpp>

namespace NetworKit {

/**
 * @ingroup matching
 * Path growing matching algorithm as described by
 * Hougardy and Drake, http://dx.doi.org/10.1016/S0020-0190(02)00393-9
 * Computes an approximate maximum weight matching with guarantee 1/2.
 * (Note that better algorithms in terms of approximation quality exist.)
 */
class PathGrowingMatcher final: public Matcher {
    void checkInput() const;

public:
    /**
     * @param[in] G Undirected graph with no self-loops for which matching is computed.
     */
    PathGrowingMatcher(const Graph& G);

    /**
     * @param[in] G Undirected graph with no self-loops for which matching is computed.
     */
    PathGrowingMatcher(const Graph& G, const std::vector<double>& edgeScores);

    /**
     * Runs path growing algorithm to compute approximate maximum weight matching
     * for graph @a G.
     * @return Matching (at least half as heavy as maximum weight matching).
     */
    void run() override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_MATCHING_PATH_GROWING_MATCHER_HPP_
