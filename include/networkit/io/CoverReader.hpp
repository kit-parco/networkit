
#ifndef NETWORKIT_IO_COVER_READER_HPP_
#define NETWORKIT_IO_COVER_READER_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class CoverReader {
public:
    /**
     * Read a cover from a file. File format: each line contains the node ids of one subset.
     *
     * @param[in]  path  The path to the input file
     * @param[in]  G     The graph for which the cover shall be read
     * @return The cover instance
     */
    virtual Cover read(const std::string &path, Graph &G);
};

} /* namespace NetworKit */

#endif // NETWORKIT_IO_COVER_READER_HPP_
