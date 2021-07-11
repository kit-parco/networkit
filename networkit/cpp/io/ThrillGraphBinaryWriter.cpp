// no-networkit-format
/*
 * ThrillGraphBinaryWriter.cpp
 *
 * @author Michael Hamann
 */

#include <fstream>

#include <networkit/io/ThrillGraphBinaryWriter.hpp>

namespace NetworKit {

void ThrillGraphBinaryWriter::write(const Graph &G, const std::string &path) {
    if (G.upperNodeIdBound() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("Thrill binary graphs only support graphs with up to 2^32-1 nodes.");
    }

    std::ofstream out_stream(path, std::ios::trunc | std::ios::binary);

    std::vector<uint32_t> neighbors;

    for (node u = 0; u < G.upperNodeIdBound(); ++u) {
        neighbors.clear();

        if (G.hasNode(u)) {
            G.forEdgesOf(u, [&](node v) {
                if (u <= v) {
                    neighbors.push_back(v);
                }
            });
        }

        size_t deg = neighbors.size();

        // Write variable length int for the degree
        if(!deg) {
            out_stream << uint8_t(0);
        }

        while(deg) {
            size_t u = deg & 0x7F;
            deg >>= 7;
            out_stream << uint8_t(u | (deg ? 0x80 : 0));
        }

        for (uint32_t v : neighbors) {
            // write neighbor as little endian
            for (size_t i = 0; i < sizeof(uint32_t); ++i) {
                out_stream << uint8_t(v);
                v >>= 8;
            }
        }
    }

    out_stream.close();
}
} // namespace NetworKit
