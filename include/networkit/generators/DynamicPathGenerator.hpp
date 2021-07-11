// no-networkit-format
/*
 * DynamicPathGenerator.hpp
 *
 *  Created on: 14.01.2014
 *      Author: cls
 */

#ifndef NETWORKIT_GENERATORS_DYNAMIC_PATH_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_DYNAMIC_PATH_GENERATOR_HPP_

#include <networkit/generators/DynamicGraphGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Example dynamic graph generator: Generates a dynamically growing path.
 */
class DynamicPathGenerator final : public DynamicGraphGenerator {
public:

    std::vector<GraphEvent> generate(count nSteps) override;

};

} /* namespace NetworKit */

#endif // NETWORKIT_GENERATORS_DYNAMIC_PATH_GENERATOR_HPP_
