
#ifndef NETWORKIT_IO_EDGE_LIST_COVER_READER_HPP_
#define NETWORKIT_IO_EDGE_LIST_COVER_READER_HPP_

#include <networkit/io/CoverReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class EdgeListCoverReader final : public CoverReader {
public:
    /**
     * Constructs the EdgeListCoverReader class with @a firstNode as the index of the first node in
     * the file.
     * @param[in] firstNode Index of the first node in the file.
     */
    EdgeListCoverReader(node firstNode = 1);

    /**
     * Read a cover from a file. File format: each line contains the node ids of one subset.
     *
     * @param[in]  path  The path to the input file
     * @param[in]  G     The graph for which the cover shall be read
     * @return The cover instance
     */
    Cover read(const std::string &path, Graph &G) override;

private:
    node firstNode;
};
} // namespace NetworKit

#endif // NETWORKIT_IO_EDGE_LIST_COVER_READER_HPP_
