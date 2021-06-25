/*
 * NetworkitBinaryReader.hpp
 *
 *      Author: Charmaine Ndolo <charmaine.ndolo@hu-berlin.de>
 */

// networkit-format

#ifndef NETWORKIT_IO_NETWORKIT_BINARY_READER_HPP_
#define NETWORKIT_IO_NETWORKIT_BINARY_READER_HPP_

#include <cstring>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphReader.hpp>
#include <networkit/io/MemoryMappedFile.hpp>

namespace NetworKit {

/**
 * @ingroup io
 *
 * Reads a graph written in the custom Networkit format documented in cpp/io/NetworkitGraph.md
 */

class NetworkitBinaryReader final : public GraphReader {

public:
    NetworkitBinaryReader(){};

    Graph read(const std::string &path) override;
    Graph readState(const std::string &data);

private:
    count nodes;
    count chunks;
    bool directed;
    bool weighted;
    bool indexed;

    template <class T>
    Graph readData(T &source);

    template <class T>
    const char *getIterator(T &source);

    template <>
    const char *getIterator(MemoryMappedFile &source) {
        return source.cbegin();
    };

    template <>
    const char *getIterator(std::string &source) {
        return source.data();
    };
};
} // namespace NetworKit

#endif // NETWORKIT_IO_NETWORKIT_BINARY_READER_HPP_
