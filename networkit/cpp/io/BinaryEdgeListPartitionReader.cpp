// no-networkit-format
#include <fstream>

#include <networkit/io/BinaryEdgeListPartitionReader.hpp>

namespace NetworKit {

BinaryEdgeListPartitionReader::BinaryEdgeListPartitionReader(node firstNode, uint8_t width) : firstNode(firstNode), width(width) {
    if (width != 4 && width != 8) {
        throw std::runtime_error("Error: width must be 4 or 8");
    }
}

Partition BinaryEdgeListPartitionReader::read(const std::string& path) {
    return read(std::vector<std::string>(1, path));
}

Partition BinaryEdgeListPartitionReader::read(const std::vector<std::string>& paths) {
    Partition zeta(0);

    if (!paths.empty()) {
        std::ifstream is;

        count file_index = 0;

        auto next_input = [&]() {
            is.close();
            if (file_index < paths.size()) {
                is.open(paths[file_index++], std::ios_base::in | std::ios_base::binary);
                is.exceptions(std::ifstream::failbit | std::ifstream::badbit);

                if (!is) {
                    throw std::runtime_error("Error: partition file couldn't be opened");
                }
            }
        };

        auto read_little_endian = [](std::ifstream &is, uint8_t width) -> index {
            index result = 0;
            for (uint8_t i = 0; i < width; ++i) {
                index t = is.get();
                result |= (t << (i * 8));
            }

            return result;
        };

        next_input();


        index read_values = 0;
        while (is.good() && is.is_open()) {
            index u = read_little_endian(is, width);
            index p = read_little_endian(is, width);

            if (u < firstNode) {
                throw std::runtime_error("Error: node smaller than the given firstNode found!");
            }

            u -= firstNode;

            if (p != none && p >= zeta.upperBound()) {
                zeta.setUpperBound(p + 1);
            }

            while (u >= zeta.numberOfElements()) {
                zeta.extend();
            }

            zeta[u] = p;

            if (is.is_open() && (is.peek() == std::char_traits<char>::eof() || !is.good())) {
                next_input();
            }

            ++read_values;
        }

        if (read_values < zeta.numberOfElements()) {
            throw std::runtime_error("Error, read less values than there are elements in the partition.");
        } else if (read_values > zeta.numberOfElements()) {
            throw std::runtime_error("Error, read more values than there are elements in the partition.");
        }
    }

    return zeta;
}

} // namespace NetworKit
