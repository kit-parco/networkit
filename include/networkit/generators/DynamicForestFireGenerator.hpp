// no-networkit-format
/*
 * DynamicForestFireGenerator.hpp
 *
 *  Created on: 17.01.2014
 *      Author: cls
 */

#ifndef NETWORKIT_GENERATORS_DYNAMIC_FOREST_FIRE_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_DYNAMIC_FOREST_FIRE_GENERATOR_HPP_

#include <networkit/generators/DynamicGraphGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * The Forest Fire generative model produces dynamic graphs with the following properties:
 * - heavy tailed degree distribution
 * - communities
 * - densification power law
 * - shrinking diameter
 *
 * see Leskovec, Kleinberg, Faloutsos: Graphs over Tim: Densification Laws,
 *   Shringking Diameters and Possible Explanations
 */
class DynamicForestFireGenerator final : public DynamicGraphGenerator {

public:

    /**
     * @param  p        "forward burning probability", controls the amount of forward connections burned by each new node
     * @param   directed  whether the generated graph is directed or undirected
     * @param   r         scales down the burning probability for backward connections (ignored for undirected graphs)
     */
    DynamicForestFireGenerator(double p, bool directed, double r = 1.0);

    /**
     * Generate event stream.
     *
     * @param[in]	nSteps	number of time steps in the event stream
     */
     std::vector<GraphEvent> generate(count nSteps) override;

private:

     double p;
     bool directed;
     double r;
     bool firstCall;
     Graph G;

};

} /* namespace NetworKit */

#endif // NETWORKIT_GENERATORS_DYNAMIC_FOREST_FIRE_GENERATOR_HPP_
