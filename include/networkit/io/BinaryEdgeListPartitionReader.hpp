
#ifndef NETWORKIT_IO_BINARY_EDGE_LIST_PARTITION_READER_HPP_
#define NETWORKIT_IO_BINARY_EDGE_LIST_PARTITION_READER_HPP_

#include <string>

#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * Reads a partition file that contains a binary list of pairs (node, partition(node)).
 * It is assumed that all integers are unsigned.
 *
 * @ingroup io
 */
class BinaryEdgeListPartitionReader final {

public:
    /**
     * Constructs the BinaryEdgeListPartitionReader class with @a firstNode as the id of the node
     * that shall become node 0 and @a width the width of the integers read.
     * @param[in] firstNode Index of the first node in the file.
     * @param[in] width The width of the read integers (supported values: 4, 8)
     */
    BinaryEdgeListPartitionReader(node firstNode = 0, uint8_t width = 4);

    /**
     * Read a partition from a file. File format:
     * binary list of pairs (node, partition(node))
     * Partition ids must be in the range [0, max(uint64_t)).
     *
     * @param[in] path Path to file or to several files (which are read in order).
     * @return The partition contained in the file at @a path.
     */
    Partition read(const std::string &path);

    Partition read(const std::vector<std::string> &paths);

    node firstNode;
    uint8_t width;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_BINARY_EDGE_LIST_PARTITION_READER_HPP_
