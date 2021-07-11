// no-networkit-format
/*
 * Matcher.hpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_MATCHING_MATCHER_HPP_
#define NETWORKIT_MATCHING_MATCHER_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/matching/Matching.hpp>

namespace NetworKit {

/**
 * @ingroup matching
 * Abstract base class for matching algorithms.
 */
class Matcher : public Algorithm {
protected:
    const Graph* G;
    Matching M;
    bool edgeScoresAsWeights;  //<! if true, algorithm should use edge scores instead of edge weights
    const std::vector<double> edgeScores;   //<! optional edge scores to be used instead of edge weight

public:
    /**
     * Constructor.
     * @param[in] G Graph for which matching is to be computed.
     */
    Matcher(const Graph& G);

    /**
     * Constructor.
     * @param[in] G Graph for which matching is to be computed.
     * @param edgeScores 	(optional) to be used instead of weights
     */
    Matcher(const Graph& G, const std::vector<double>& edgeScores);

    /** Default destructor */
    ~Matcher() override = default;

    /**
     * Run the matching algorithm on the stored graph and return a matching.
     * @return A matching of the stored graph.
     */
    void run() override = 0;

    Matching getMatching() const;
};

} /* namespace NetworKit */
#endif // NETWORKIT_MATCHING_MATCHER_HPP_
