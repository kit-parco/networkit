// no-networkit-format
/*
 * DynamicGraphGenerator.hpp
 *
 *  Created on: 14.01.2014
 *      Author: cls
 */

#ifndef NETWORKIT_GENERATORS_DYNAMIC_GRAPH_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_DYNAMIC_GRAPH_GENERATOR_HPP_

#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Abstract base class for a dynamic graph generator (in the new dynamic architecture).
 * The generator produces a stream of events.
 */
class DynamicGraphGenerator {

public:
    /** Default destructor */
    virtual ~DynamicGraphGenerator() = default;

    /**
     * Generate event stream.
     *
     * @param[in]	nSteps	number of time steps in the event stream
     */
    virtual std::vector<GraphEvent> generate(count nSteps) = 0;

protected:

    Graph G; // the graph instance
};

} /* namespace NetworKit */

#endif // NETWORKIT_GENERATORS_DYNAMIC_GRAPH_GENERATOR_HPP_
