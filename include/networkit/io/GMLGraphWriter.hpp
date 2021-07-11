/*
 * GMLGraphWriter.hpp
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_GML_GRAPH_WRITER_HPP_
#define NETWORKIT_IO_GML_GRAPH_WRITER_HPP_

#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class GMLGraphWriter final : public GraphWriter {
public:
    /** Default constructor */
    GMLGraphWriter() = default;

    /**
     * Write a graph @a G and its coordinates to a GML file.
     *
     * @param[in]  G     Graph of type NetworKit with 2D coordinates.
     * @param[in]  path  Path to file.
     */
    void write(const Graph &G, const std::string &path) override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_GML_GRAPH_WRITER_HPP_
