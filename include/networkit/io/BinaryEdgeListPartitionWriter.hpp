/**
 * BinaryEdgeListPartitionWriter.hpp
 *
 * @author Michael Hamann
 */

#ifndef NETWORKIT_IO_BINARY_EDGE_LIST_PARTITION_WRITER_HPP_
#define NETWORKIT_IO_BINARY_EDGE_LIST_PARTITION_WRITER_HPP_

#include <string>

#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Writes a partition to a file that contains a binary list of pairs (node, partition(node)).
 * All integers are unsigned integers.
 */
class BinaryEdgeListPartitionWriter final {
public:
    /**
     * Constructs the BinaryEdgeListPartitionWriter class with @a firstNode as the id of the
     * first node (i.e., this value is added to all node ids before writing) using unsigned
     * integers with the given @a width.
     *
     * @param[in] firstNode The id of node 0, all other node ids are adjusted accordingly.
     * @param[in] width The width of the written integers (supported values: 4, 8).
     */
    BinaryEdgeListPartitionWriter(node firstNode = 0, uint8_t width = 4);

    /**
     * Write the given partition @a zeta to the given @a path.
     *
     * @param[in] zeta The partition to write.
     * @param[in] path The path to write to.
     */
    void write(Partition &zeta, const std::string &path) const;

private:
    node firstNode;
    uint8_t width;
};
} /* namespace NetworKit */

#endif // NETWORKIT_IO_BINARY_EDGE_LIST_PARTITION_WRITER_HPP_
