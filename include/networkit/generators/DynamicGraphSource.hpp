// no-networkit-format
/*
 * DynamicGraphSpurce.hpp
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#ifndef NETWORKIT_GENERATORS_DYNAMIC_GRAPH_SOURCE_HPP_
#define NETWORKIT_GENERATORS_DYNAMIC_GRAPH_SOURCE_HPP_

#include <functional>

#include <networkit/dynamics/GraphEventProxy.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 */
class DynamicGraphSource {

protected:

    GraphEventProxy* Gproxy; //!< receives events produced by the generator and forwards them
    Graph* G;
    bool graphSet; //!< true if newGraph has been called and graph and proxy instances are properly set
    bool graphInitialized; //!< true if initializeGraph has been called and graph has been properly initialized

public:

    /** Default constructor */
     DynamicGraphSource();

    // DynamicGraphGenerator(GraphEventProxy& proxy);

    virtual ~DynamicGraphSource() = default;

    /**
     * After constructing a DynamicGraphGenerator, call this to set a new
     * a Graph and GraphEventProxy instance and get access to them.
     */
    GraphEventProxy* newGraph();

    /**
     * The generator may expect the graph to be in a certain initial state. Call this method first.
     */
    virtual void initializeGraph() = 0;

    /**
     * Perform one generative step - as defined by the implementation.
     */
    virtual void generate() = 0;

    /*
     * Continue generating while function does not return false.
     * @param[in] cont generator continues when this function returns true
     */
    virtual void generateWhile(std::function<bool(void)> cont);

    /**
     * Continue generating until the number of nodes reaches this upper limit.
     * @param[in] n  umber of nodes
     */
    virtual void generateNodes(count n);

    /**
     * Continue generating until the number of edges reaches this upper limit.
     * @param[in] m  umber of edges
     */
    virtual void generateEdges(count m);

    /**
     * Continue generating until the number of time steps reaches this upper limit.
     */
    virtual void generateTimeSteps(count t);
};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_DYNAMIC_GRAPH_SOURCE_HPP_
