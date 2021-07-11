// no-networkit-format
#include <algorithm>
#include <fstream>
#include <stdexcept>

#include <networkit/graph/GraphBuilder.hpp>
#include <networkit/io/ThrillGraphBinaryReader.hpp>

namespace NetworKit {

ThrillGraphBinaryReader::ThrillGraphBinaryReader(count n) : n(n) {}

namespace {
    uint32_t get_uint32(std::ifstream& is) {
        uint32_t result = 0;

        for (size_t i = 0; i < sizeof(uint32_t); ++i) {
            uint32_t u = is.get();
            result |= (u << (i * 8));
        }

        return result;
    }

    uint64_t get_variant(std::ifstream& is) {
        size_t v = 0;
        for (size_t shift_width = 0; shift_width < 64; shift_width += 7) {
            uint64_t u = is.get();

            // The last value is just a single bit (the 64th bit) - throw if there is
            // more than this single bit in the input.
            if (shift_width == 63 && (u & 0xFE)) {
                throw std::overflow_error("Overflow during variant64 decoding.");
            }

            // Read exactly 7 bits from the input
            v |= (u & 0x7F) << shift_width;

            // If the 8th bit is not set, we reached the end of the variant
            if (!(u & 0x80)) break;
        }

        return v;
    }
}


Graph ThrillGraphBinaryReader::read(const std::string &path) {
    return read(std::vector<std::string>(1, path));
}

Graph ThrillGraphBinaryReader::read(const std::vector<std::string> &paths) {
    GraphBuilder gb(n);

    if (!paths.empty()) {
        std::ifstream is;
        count file_index = 0;

        auto next_input = [&]() {
            is.close();
            if (file_index < paths.size()) {
                is.open(paths[file_index++], std::ios::binary);
                // Throw an exception when no byte could be read at any point
                is.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            }
        };

        next_input();

        node max_id = 0;

        for (node u = 0; is.good() && is.is_open(); ++u) {
            // Add node if it does not exist yet, only one may be missing
            if (u >= gb.upperNodeIdBound()) {
                gb.addNode();
            }

            for (count deg = get_variant(is); deg > 0 && is.good(); --deg) {
                uint32_t v = get_uint32(is);

                max_id = std::max<node>(max_id, v);

                gb.addHalfEdge(u, v);
            }

            if (is.is_open() && (is.peek() == std::char_traits<char>::eof() || !is.good())) {
                next_input();
            }
        }

        if (max_id >= gb.upperNodeIdBound()) {
            throw std::runtime_error("Maximum read node id larger than number of nodes read.");
        }
    }

    return gb.toGraph(true, true);
}

} // namespace NetworKit
