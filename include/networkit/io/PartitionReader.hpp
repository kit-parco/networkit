/*
 * PartitionReader.hpp
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_PARTITION_READER_HPP_
#define NETWORKIT_IO_PARTITION_READER_HPP_

#include <fstream>

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class PartitionReader final {

public:
    /**
     * Read a clustering from a file. File format:
     *     line n contains cluster id of node (n - 1)
     *
     * @param[in]  path  Path to file.
     */
    Partition read(const std::string &path);
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_PARTITION_READER_HPP_
