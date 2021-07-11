/*
 * StronglyConnectedComponents.hpp
 *
 *  Created on: 01.06.2014
 *      Authors: Klara Reichard <klara.reichard@gmail.com>
 *               Marvin Ritter <marvin.ritter@gmail.com>
 *               Obada Mahdi <omahdi@gmail.com>
 *               Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_COMPONENTS_STRONGLY_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_STRONGLY_CONNECTED_COMPONENTS_HPP_

#include <networkit/components/ComponentDecomposition.hpp>

namespace NetworKit {

/**
 * @ingroup components
 */
class StronglyConnectedComponents final : public ComponentDecomposition {

public:
    /**
     * Determines the strongly connected components of a directed graph using Tarjan's algorithm.
     *
     * @param G Graph A directed graph.
     */
    StronglyConnectedComponents(const Graph &G);

    /**
     * Runs the algorithm.
     */
    void run() override;

    /**
     * Constructs a new graph that contains only the nodes inside the largest
     * strongly connected component.
     * @param G            The input graph.
     * @param compactGraph If true, the node ids of the output graph will be compacted
     * (i.e. re-numbered from 0 to n-1). If false, the node ids will not be changed.
     */
    static Graph extractLargestStronglyConnectedComponent(const Graph &G,
                                                          bool compactGraph = false);
};

} // namespace NetworKit

#endif // NETWORKIT_COMPONENTS_STRONGLY_CONNECTED_COMPONENTS_HPP_
