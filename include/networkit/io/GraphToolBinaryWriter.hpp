/*
 * GraphToolBinaryWriter.hpp
 *
 *  Created on: 02.12.14
 *      Author: Maximilian Vogel
 */

#ifndef NETWORKIT_IO_GRAPH_TOOL_BINARY_WRITER_HPP_
#define NETWORKIT_IO_GRAPH_TOOL_BINARY_WRITER_HPP_

#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**
 * Writes graphs to files in the binary format defined by graph-tool[1].
 * [1]: http://graph-tool.skewed.de/static/doc/gt_format.html
 */
class GraphToolBinaryWriter final : public GraphWriter {

public:
    GraphToolBinaryWriter(bool littleEndianness = true);

    /**
     * Given the path of an input file, read the graph contained.
     *
     * @param[in]  path  input file path
     */
    void write(const Graph &G, const std::string &path) override;

private:
    bool littleEndianness;

    void writeAdjacencies(std::ofstream &file, const Graph &G);

    uint8_t getAdjacencyWidth(uint64_t n);

    void writeComment(std::ofstream &file);

    void writeHeader(std::ofstream &file);

    template <typename Type>
    void writeType(std::ofstream &file, int width, Type val) {
        // DEBUG("writing ",val, "with width ", width, " to file");
        uint8_t *bytes = new uint8_t[width];
        if (this->littleEndianness) {
            for (int i = 0; i < width; ++i) {
                bytes[i] = (val >> (i * 8)) & 0xFF;
            }
        } else {
            for (int i = 0; i < width; ++i) {
                bytes[i] = (val >> ((width - 1 - i) * 8)) & 0xFF;
            }
        }
        file.write((char *)bytes, width);
        delete[] bytes;
    }
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_GRAPH_TOOL_BINARY_WRITER_HPP_
