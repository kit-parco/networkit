// no-networkit-format
/*
 * RandomSpanningForest.hpp
 *
 *  Created on: 06.09.2015
 *      Author: Henning
 */

#ifndef NETWORKIT_COMPONENTS_RANDOM_SPANNING_FOREST_HPP_
#define NETWORKIT_COMPONENTS_RANDOM_SPANNING_FOREST_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/graph/SpanningForest.hpp>

namespace NetworKit {

/**
 * Creates a random spanning tree for each connected component.
 * Time complexity: cover time of G.
 * @ingroup graph
 */
class RandomSpanningForest final : public SpanningForest {
public:
    RandomSpanningForest(const Graph& G);

    ~RandomSpanningForest() override = default;

    /**
     * Computes for each component a random spanning tree.
     * Uses simple random-walk based algorithm.
     * Time complexity: cover time of G.
     */
    void run() override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMPONENTS_RANDOM_SPANNING_FOREST_HPP_
