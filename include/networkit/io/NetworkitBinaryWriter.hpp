/*
 * NetworkitBinaryWriter.hpp
 *
 *      Author: Charmaine Ndolo <charmaine.ndolo@hu-berlin.de>
 */

// networkit-format

#ifndef NETWORKIT_IO_NETWORKIT_BINARY_WRITER_HPP_
#define NETWORKIT_IO_NETWORKIT_BINARY_WRITER_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

enum class NetworkitBinaryWeights {
    none,
    unsignedFormat,
    signedFormat,
    doubleFormat,
    floatFormat,
    autoDetect
};

/**
 * @ingroup io
 *
 * Writes a graph in the custom Networkit format documented in cpp/io/NetworkitGraph.md
 */
class NetworkitBinaryWriter final : public GraphWriter {

public:
    NetworkitBinaryWriter(uint64_t chunks = 32,
                          NetworkitBinaryWeights weightsType = NetworkitBinaryWeights::autoDetect,
                          bool preserveEdgeIndex = false);

    void write(const Graph &G, const std::string &path) override;
    std::string writeState(const Graph &G);

private:
    count chunks;
    NetworkitBinaryWeights weightsType;
    bool preserveEdgeIndex;

    template <class T>
    void writeData(T &outStream, const Graph &G);
};

} // namespace NetworKit

#endif // NETWORKIT_IO_NETWORKIT_BINARY_WRITER_HPP_
