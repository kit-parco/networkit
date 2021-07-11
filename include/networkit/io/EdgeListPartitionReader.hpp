
#ifndef NETWORKIT_IO_EDGE_LIST_PARTITION_READER_HPP_
#define NETWORKIT_IO_EDGE_LIST_PARTITION_READER_HPP_

#include <fstream>

#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class EdgeListPartitionReader final {

public:
    /**
     * Constructs the EdgeListPartitionReader class with @a firstNode as the index of the first node
     * in the file.
     * @param[in]  firstNode  Index of the first node in the file.
     * @param[in]  sepChar    The separator between two elements
     */
    EdgeListPartitionReader(node firstNode = 1, char sepChar = '\t');

    /**
     * Read a clustering from a file. File format:
     *     A list of the nodes and their membership (memberships are labelled by integer numbers
     * >=1).
     *
     * @param[in]  path  Path to file.
     * @return The clustering contained in the file at @a path.
     */
    Partition read(const std::string &path);

private:
    node firstNode;
    char sepChar;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_EDGE_LIST_PARTITION_READER_HPP_
