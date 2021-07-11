// no-networkit-format
#include <fstream>

#include <networkit/io/BinaryPartitionReader.hpp>

namespace NetworKit {

BinaryPartitionReader::BinaryPartitionReader(uint8_t width) : width(width) {
    if (width != 4 && width != 8) {
        throw std::runtime_error("Only 4 and 8 are supported widths");
    }
}

Partition BinaryPartitionReader::read(const std::string& path) {
        std::ifstream is(path, std::ios_base::in | std::ios_base::binary);

        if (!is) {
            throw std::runtime_error("Error: partition file couldn't be opened");
        }

        is.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        is.seekg(0, std::ios_base::end);
        if ((is.tellg() % width) != 0) {
            throw std::runtime_error("Error: length of partition file must be a multiple of the width.");
        }

        count length = is.tellg() / width;
        is.seekg(0);

        Partition zeta(length);

        for (index u = 0; u < length; ++u) {
            index p = 0;

            for (size_t i = 0; i < width; ++i) {
                uint64_t t = is.get();
                p |= (t << (i * 8));
            }

            if (p != none && p >= zeta.upperBound()) {
                zeta.setUpperBound(p + 1);
            }

            zeta[u] = p;
        }

        return zeta;
}

} // namespace NetworKit
