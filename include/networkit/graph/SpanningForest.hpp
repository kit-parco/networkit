// no-networkit-format
/*
 * SpanningForest.hpp
 *
 *  Created on: 06.09.2015
 *      Author: Henning
 */

#ifndef NETWORKIT_GRAPH_SPANNING_FOREST_HPP_
#define NETWORKIT_GRAPH_SPANNING_FOREST_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Base class for spanning forest/tree algorithms.
 */
class SpanningForest {
protected:
    const Graph* G;
    Graph forest;

public:
    SpanningForest(const Graph& G) : G(&G) {}
    virtual ~SpanningForest() = default;

    virtual void run();

    /**
     * @return Forest computed by run method.
     * Note: So far no explicit check if run method has been invoked before.
     */
    Graph getForest() {
        return forest;
    }
};

} /* namespace NetworKit */
#endif // NETWORKIT_GRAPH_SPANNING_FOREST_HPP_
