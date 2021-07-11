// no-networkit-format
#include <fstream>

#include <networkit/io/BinaryEdgeListPartitionWriter.hpp>

namespace NetworKit {

BinaryEdgeListPartitionWriter::BinaryEdgeListPartitionWriter(node firstNode, uint8_t width) : firstNode(firstNode), width(width) {
    if (width != 4 && width != 8) {
        throw std::runtime_error("Width must be 4 or 8");
    }
}

void BinaryEdgeListPartitionWriter::write( Partition &zeta, const std::string &path ) const {
    auto write_little_endian = [](std::ofstream &os, index x, uint8_t width) {
        for (uint8_t w = 0; w < width; ++w) {
            os.put(uint8_t(x));
            x >>= 8;
        }
    };

    if (width == 4 && zeta.upperBound() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("Error, the upper bound of the given partition cannot be represented by an unsigned int of width 4. Please use a width of 8.");
    }

    std::ofstream os(path, std::ios::trunc | std::ios::binary);

    os.exceptions(std::ofstream::badbit | std::ofstream::failbit);

    zeta.forEntries([&](index u, index p) {
        write_little_endian(os, u + firstNode, width);
        write_little_endian(os, p, width);
    });
}

} // namespace NetworKit
