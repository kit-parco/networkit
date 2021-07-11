// no-networkit-format
/*
 * ErdosRenyiGenerator.hpp
 *
 *  Created on: 21.01.2014
 *      Author: Henning
 */

#ifndef NETWORKIT_GENERATORS_ERDOS_RENYI_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_ERDOS_RENYI_GENERATOR_HPP_

#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {
/**
 * @ingroup generators
 * Creates G(n, p) graphs.
 *
 * This class is a wrapper to @ref ErdosRenyiEnumerator.
 */
class ErdosRenyiGenerator final : public StaticGraphGenerator {
public:
    /**
     * Creates random graphs in the G(n,p) model.
     * The generation follows a parallelized version of Vladimir Batagelj
     * and Ulrik Brandes: "Efficient generation of large random networks",
     * Phys Rev E 71, 036113 (2005) with a runtime of O((n+m) / P) (WHP),
     * where n, m, and P are the numbers of nodes, edges and parallel threads
     * respectively.
     *
     * The generator can also be used to generate empty or complete graphs
     * efficiently by setting (p = 0.0 or p = 1.0). Observe that this is
     * implemented as a special case and requires p to be exactly zero or
     * one.
     *
     * @warning For compatibility reasons, the generator does not produce
     * self-loops by default.
     *
     * @param nNodes Number of nodes n in the graph.
     * @param prob Probability of existence for each edge p.
     * @param directed  generates a directed graph
     * @param self_loops Controls whether a directed graph may contain self_loops
     *                   (undirected graphs never have them)
     */
    ErdosRenyiGenerator(count nNodes, double prob, bool directed = false, bool self_loops = false);

    Graph generate() override;

private:
    node nNodes;
    double prob;
    bool directed;
    bool self_loops;
};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_ERDOS_RENYI_GENERATOR_HPP_
