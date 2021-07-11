/*
 * ThrillGraphBinaryReader.hpp
 *
 * @author Michael Hamann
 */

#ifndef NETWORKIT_IO_THRILL_GRAPH_BINARY_READER_HPP_
#define NETWORKIT_IO_THRILL_GRAPH_BINARY_READER_HPP_

#include <string>
#include <vector>

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * Reads a graph format consisting of a serialized DIA of vector<uint32_t> from thrill.
 */
class ThrillGraphBinaryReader final : public GraphReader {

public:
    /**
     * When the number of nodes is given, reading the graph is more efficient.
     *
     * @param[in]  n  The number of nodes
     */
    ThrillGraphBinaryReader(count n = 0);

    /**
     * Given the path of an input file, read the graph contained.
     *
     * @param[in]  path  input file path
     */
    Graph read(const std::string &path) override;

    Graph read(const std::vector<std::string> &path);

private:
    const count n;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_THRILL_GRAPH_BINARY_READER_HPP_
