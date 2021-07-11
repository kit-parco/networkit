/*
 * GraphWriter.hpp
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_GRAPH_WRITER_HPP_
#define NETWORKIT_IO_GRAPH_WRITER_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Abstract base class for graph writers.
 */
class GraphWriter {
public:
    virtual ~GraphWriter() = default;

    virtual void write(const Graph &G, const std::string &path) = 0;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_GRAPH_WRITER_HPP_
