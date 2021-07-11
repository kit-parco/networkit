#ifndef NETWORKIT_SCD_TWO_PHASE_L_HPP_
#define NETWORKIT_SCD_TWO_PHASE_L_HPP_

#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {

/**
 * The two-phase local community detection algorithm optimizing the L-measure.
 *
 * This is an implementation of the algorithm proposed in:
 *
 * Chen, J., Zaïane, O., & Goebel, R. (2009).
 * Local Community Identification in Social Networks.
 * In 2009 International Conference on Advances in Social Network Analysis and Mining (pp. 237–242).
 * https://doi.org/10.1109/ASONAM.2009.14
 */
class TwoPhaseL : public SelectiveCommunityDetector {

public:
    /**
     * Construct the algorithm class.
     *
     * @param G The graph on which communities shall be found.
     */
    TwoPhaseL(const Graph &g);

    /**
     * @param[in]	seeds	seed nodes
     *
     * @param[out]		community as a set of nodes
     */
    std::set<node> expandOneCommunity(const std::set<node> &s) override;

    // inherit method from parent class.
    using SelectiveCommunityDetector::expandOneCommunity;
};

} /* namespace NetworKit */
#endif // NETWORKIT_SCD_TWO_PHASE_L_HPP_
