/*
 * BinaryPartitionReader.hpp
 *
 *  Created on: 12.04.2017
 *      Author: Michael Hamann
 */

#ifndef NETWORKIT_IO_BINARY_PARTITION_READER_HPP_
#define NETWORKIT_IO_BINARY_PARTITION_READER_HPP_

#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class BinaryPartitionReader final {

public:
    /**
     * Construct a binary partition reader.
     *
     * @param[in]  width  The integer width. Supported values: 4 and 8.
     */
    BinaryPartitionReader(uint8_t width = 4);

    /**
     * Read a partition from a file. File format:
     *     list of (unsigned) integer partition ids, one for every node
     *
     * @param[in]  path  Path to file.
     */
    Partition read(const std::string &path);

private:
    uint8_t width;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_BINARY_PARTITION_READER_HPP_
