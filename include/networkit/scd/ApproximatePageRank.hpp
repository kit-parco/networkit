/*
 * ApproximatePageRank.hpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#ifndef NETWORKIT_SCD_APPROXIMATE_PAGE_RANK_HPP_
#define NETWORKIT_SCD_APPROXIMATE_PAGE_RANK_HPP_

#include <set>
#include <unordered_map>
#include <vector>

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Computes an approximate PageRank vector from a given seed.
 */
class ApproximatePageRank final {
    const Graph *g;
    double alpha;
    double eps;

    std::unordered_map<node, std::pair<double, double>> prRes;

public:
    /**
     * @param g Graph for which an APR is computed.
     * @param alpha Loop probability of random walk.
     * @param epsilon Error tolerance.
     */
    ApproximatePageRank(const Graph &g, double alpha, double epsilon = 1e-12);

    /**
     * @return Approximate PageRank vector from @a seeds with parameters
     *         specified in the constructor.
     */
    std::vector<std::pair<node, double>> run(const std::set<node> &seeds);

    /**
     * @return Approximate PageRank vector from @a seed with parameters
     *         specified in the constructor.
     */
    std::vector<std::pair<node, double>> run(node seed);
};

} /* namespace NetworKit */
#endif // NETWORKIT_SCD_APPROXIMATE_PAGE_RANK_HPP_
