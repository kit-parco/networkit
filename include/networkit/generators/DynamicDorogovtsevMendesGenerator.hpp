// no-networkit-format
/*
 * DynamicDorogovtsevMendesGenerator.hpp
 *
 *  Created on: 03.02.2014
 *      Author: cls
 */

#ifndef NETWORKIT_GENERATORS_DYNAMIC_DOROGOVTSEV_MENDES_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_DYNAMIC_DOROGOVTSEV_MENDES_GENERATOR_HPP_

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/DynamicGraphGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 */
class DynamicDorogovtsevMendesGenerator final: public DynamicGraphGenerator {

public:

    DynamicDorogovtsevMendesGenerator();

    std::vector<GraphEvent> generate(count nSteps) override;

private:

    std::vector<std::pair<node, node>> edges;
    bool initial;
    node u; // current node

};

} /* namespace NetworKit */

#endif // NETWORKIT_GENERATORS_DYNAMIC_DOROGOVTSEV_MENDES_GENERATOR_HPP_
