/*
 * GlobalThresholdFilter.hpp
 *
 *  Created on: 23.07.2014
 *      Author: Gerd Lindner
 */

#ifndef NETWORKIT_SPARSIFICATION_GLOBAL_THRESHOLD_FILTER_HPP_
#define NETWORKIT_SPARSIFICATION_GLOBAL_THRESHOLD_FILTER_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Calculates a sparsified graph by applying a global threshold to an edge score.
 */
class GlobalThresholdFilter {
    const Graph *graph;
    const std::vector<double> &attribute;
    const double threshold;
    const bool above;

public:
    /**
     * Creates a new instance of a global threshold filter.
     * @param threshold the threshold
     * @param above if set to true, edges with a score above or equal to the threshold remain in
     * the filtered graph. If set to false, edges with a score below or equal to the threshold
     * remain in the filtered graph.
     */
    GlobalThresholdFilter(const Graph &graph, const std::vector<double> &attribute,
                          double threshold, bool above);

    Graph calculate();
};

} // namespace NetworKit
#endif // NETWORKIT_SPARSIFICATION_GLOBAL_THRESHOLD_FILTER_HPP_
