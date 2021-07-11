/**
 * BinaryPartitionWriter.hpp
 *
 * @author Michael Hamann
 */

#ifndef NETWORKIT_IO_BINARY_PARTITION_WRITER_HPP_
#define NETWORKIT_IO_BINARY_PARTITION_WRITER_HPP_

#include <string>

#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Writes a partition to a file that contains a binary list of partition ids.
 * Partition ids are unsigned integers.
 */
class BinaryPartitionWriter final {
public:
    /**
     * Constructs the BinaryPartitionWriter class using unsigned integers
     * of width @a width for the written partition ids.
     *
     * @param[in] width The width of the written integers (supported values: 4, 8, default: 4).
     */
    BinaryPartitionWriter(uint8_t width = 4);

    /**
     * Write the given partition @a zeta to the given @a path.
     *
     * @param[in] zeta The partition to write.
     * @param[in] path The path to write to.
     */
    void write(const Partition &zeta, const std::string &path) const;

private:
    uint8_t width;
};
} // namespace NetworKit

#endif // NETWORKIT_IO_BINARY_PARTITION_WRITER_HPP_
