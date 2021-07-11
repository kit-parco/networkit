/*
 * EdgeListWriter.hpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef NETWORKIT_IO_EDGE_LIST_WRITER_HPP_
#define NETWORKIT_IO_EDGE_LIST_WRITER_HPP_

#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * A writer for the edge list format. The output will contain one edge per line,
 * in the format fromNodeSEPARATORtoNode, where separator can be specified by
 * the user.
 */
class EdgeListWriter final : public GraphWriter {

public:
    EdgeListWriter() = default; // nullary constructor for Python shell

    /**
     * @param[in]  separator  character used to separate nodes in an edge line
     * @param[in]  firstNode  index of the first node in the file
     * @param[in]  bothDirections  for undirected graphs: if every edge shall be written in both
     * directions (default: false)
     */
    EdgeListWriter(char separator, node firstNode, bool bothDirections = false);

    /**
     * Write the graph to a file.
     * @param[in]  G    the graph
     * @param[in]  path  the output file path
     */
    void write(const Graph &G, const std::string &path) override;

private:
    char separator; //!< character separating nodes in an edge line
    node firstNode;
    bool bothDirections;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_EDGE_LIST_WRITER_HPP_
