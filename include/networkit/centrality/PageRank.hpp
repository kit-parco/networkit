/*
 * PageRank.h
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

#ifndef NETWORKIT_CENTRALITY_PAGE_RANK_HPP_
#define NETWORKIT_CENTRALITY_PAGE_RANK_HPP_

#include <atomic>
#include <limits>
#include <memory>

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 * Compute PageRank as node centrality measure.
 * NOTE: There is an inconsistency in the definition in Newman's book (Ch. 7) regarding
 * directed graphs; we follow the verbal description, which requires to sum over the incoming
 * edges (as opposed to outgoing ones).
 */
class PageRank final : public Centrality {

public:
    enum Norm { L1Norm, L2Norm };

    /**
     * Constructs the PageRank class for the Graph @a G
     *
     * @param[in] G Graph to be processed.
     * @param[in] damp Damping factor of the PageRank algorithm.
     * @param[in] tol Error tolerance for PageRank iteration.
     */
    PageRank(const Graph &G, double damp = 0.85, double tol = 1e-8, bool normalized = false);

    /**
     * Computes page rank on the graph passed in constructor.
     */
    void run() override;

    /**
     * Returns upper bound on the page rank: 1.0. This could be tighter by assuming e.g. a star
     * graph with n nodes.
     */
    double maximum() override;

    /**
     * Return the number of iterations performed by the algorithm.
     *
     * @return Number of iterations performed by the algorithm.
     */
    count numberOfIterations() const {
        assureFinished();
        return iterations;
    }

    // Maximum number of iterations allowed
    count maxIterations = std::numeric_limits<count>::max();

    // Norm used as stopping criterion
    Norm norm = Norm::L2Norm;

private:
    double damp;
    double tol;
    count iterations;
    bool normalized;
    std::atomic<double> max;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_PAGE_RANK_HPP_
