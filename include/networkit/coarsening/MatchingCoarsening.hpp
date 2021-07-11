// no-networkit-format
/*
 * MatchingCoarsening.hpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_COARSENING_MATCHING_COARSENING_HPP_
#define NETWORKIT_COARSENING_MATCHING_COARSENING_HPP_

#include <networkit/coarsening/GraphCoarsening.hpp>
#include <networkit/matching/Matching.hpp>

namespace NetworKit {

/**
 * @ingroup coarsening
 * Coarsens graph according to a matching.
 */
class MatchingCoarsening final : public GraphCoarsening {

public:
    MatchingCoarsening(const Graph& G, const Matching& M, bool noSelfLoops = false);

    /**
     * Contracts graph according to a matching.
     *
     * @param[in]	G	fine graph
     * @param[in]	M	matching
     * @param[in]	noSelfLoops  if true, self-loops are not produced
     *
     * @return		coarse graph
     */
    void run() override;

private:
    const Matching& M;
    bool noSelfLoops;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COARSENING_MATCHING_COARSENING_HPP_
