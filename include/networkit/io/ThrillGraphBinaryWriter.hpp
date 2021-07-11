/*
 * ThrillGraphBinaryWriter.hpp
 *
 * @author Michael Hamann
 */

#ifndef NETWORKIT_IO_THRILL_GRAPH_BINARY_WRITER_HPP_
#define NETWORKIT_IO_THRILL_GRAPH_BINARY_WRITER_HPP_

#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

class ThrillGraphBinaryWriter final : public GraphWriter {
public:
    /**
     * Write the given graph into a binary file at the given path.
     *
     * @param[in] G The graph to write.
     * @param[in] path The path where to write the graph.
     */
    void write(const Graph &G, const std::string &path) override;
};

} /* namespace NetworKit */

#endif // NETWORKIT_IO_THRILL_GRAPH_BINARY_WRITER_HPP_
