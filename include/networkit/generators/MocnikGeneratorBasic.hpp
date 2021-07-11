// no-networkit-format
/*
 * MocnikGeneratorBasic.hpp
 *
 * Created on: July 7, 2018
 * Author: Franz-Benjamin Mocnik <mail@mocnik-science.net>
 */

#ifndef NETWORKIT_GENERATORS_MOCNIK_GENERATOR_BASIC_HPP_
#define NETWORKIT_GENERATORS_MOCNIK_GENERATOR_BASIC_HPP_

#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 */
class MocnikGeneratorBasic: public StaticGraphGenerator {
    // GENERAL DATA

    /**
     * Position of each node in space.  The index of the vector is also the number of
     * the node.
     */
    std::vector<std::vector<double>> nodePositions;

    count dim;
    count n;
    double k;

public:
    /**
     * Creates random spatial graphs according to the Mocnik model.
     *
     * Please cite the following publications, in which you will find a
     * description of the model:
     *
     * Franz-Benjamin Mocnik: "The Polynomial Volume Law of Complex Networks in
     * the Context of Local and Global Optimization", Scientific Reports 8(11274)
     * 2018. doi: 10.1038/s41598-018-29131-0
     *
     * Franz-Benjamin Mocnik, Andrew Frank: "Modelling Spatial Structures",
     * Proceedings of the 12th Conference on Spatial Information Theory (COSIT),
     * 2015, pages 44-64. doi: 10.1007/978-3-319-23374-1_3
     *
     * Non-improved algorithm.
     *
     * @param dim  Dimension of the space.
     * @param n  Number of nodes in the graph.
     * @param k  Density parameter, determining the ratio of edges to nodes.
     */
    MocnikGeneratorBasic(count dim, count n, double k);

    Graph generate() override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_MOCNIK_GENERATOR_BASIC_HPP_
