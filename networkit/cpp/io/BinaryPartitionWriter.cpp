// no-networkit-format
#include <fstream>

#include <networkit/io/BinaryPartitionWriter.hpp>

namespace NetworKit {

BinaryPartitionWriter::BinaryPartitionWriter(uint8_t width) : width(width) {
    if (width != 4 && width != 8) {
        throw std::runtime_error("Only width 4 and 8 are supported");
    }
}

void BinaryPartitionWriter::write(const Partition &zeta, const std::string &path) const {
    if (width == 4 && zeta.upperBound() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("Error, the upper bound of the given partition cannot be represented by an unsigned int of width 4. Please use a width of 8.");
    }

    std::ofstream os(path, std::ios::trunc | std::ios::binary);

    os.exceptions(std::ofstream::badbit | std::ofstream::failbit);

    for (index i = 0; i < zeta.numberOfElements(); ++i) {
        index p = zeta[i];

        for (uint8_t w = 0; w < width; ++w) {
            os.put(uint8_t(p));
            p >>= 8;
        }
    }
}

} // namespace NetworKit
