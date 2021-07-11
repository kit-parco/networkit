/*
 * PartitionWriter.hpp
 *
 *  Created on: 22.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_PARTITION_WRITER_HPP_
#define NETWORKIT_IO_PARTITION_WRITER_HPP_

#include <fstream>

#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Write a clustering to a file.
 */
class PartitionWriter final {

public:
    void write(const Partition &zeta, const std::string &path) const {
        std::ofstream file{path};
        zeta.forEntries([&](node, const index c) { file << c << '\n'; });
    }
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_PARTITION_WRITER_HPP_
