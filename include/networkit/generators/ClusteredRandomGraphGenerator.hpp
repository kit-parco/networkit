/*
 * ClusteredRandomGraphGenerator.hpp
 *
 *  Created on: 28.02.2014
 *      Author: cls
 *              Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_GENERATORS_CLUSTERED_RANDOM_GRAPH_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_CLUSTERED_RANDOM_GRAPH_GENERATOR_HPP_

#include <networkit/generators/StaticGraphGenerator.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * The ClusteredRandomGraphGenerator class is used to create a clustered random graph.
 * The number of nodes and the number of edges are adjustable as well as the probabilities
 * for intra-cluster and inter-cluster edges.
 * In parallel the generated graph is not deterministic. To ensure determinism, use a single thread.
 */
class ClusteredRandomGraphGenerator final : public StaticGraphGenerator {
public:
    /**
     * Creates a clustered random graph:
     *
     * @param[in]	n	number of nodes
     * @param[in]	k	number of clusters
     * @param[in]	pIntra		intra-cluster edge probability
     * @param[in]	pInter	inter-cluster edge probability
     */
    ClusteredRandomGraphGenerator(count n, count k, double pIntra, double pInter);

    /**
     * Generates a clustered random graph with the properties given in the constructor.
     * @return The generated graph.
     */
    Graph generate() override;

    /**
     * Returns the generated ground truth communities.
     * @return The generated partition
     */
    Partition getCommunities();

private:
    const count n, k;
    const double pIntra, pInter;
    Partition zeta;
};

} /* namespace NetworKit */

#endif // NETWORKIT_GENERATORS_CLUSTERED_RANDOM_GRAPH_GENERATOR_HPP_
