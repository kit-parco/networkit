// no-networkit-format
/*
 * StaticGraphGenerator.hpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef NETWORKIT_GENERATORS_STATIC_GRAPH_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_STATIC_GRAPH_GENERATOR_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Abstract base class for static graph generators.
 */
class StaticGraphGenerator {

public:

    /** Default destructor */
    virtual ~StaticGraphGenerator() = default;

    virtual Graph generate() = 0;
};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_STATIC_GRAPH_GENERATOR_HPP_
