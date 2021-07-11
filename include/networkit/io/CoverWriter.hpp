
#ifndef NETWORKIT_IO_COVER_WRITER_HPP_
#define NETWORKIT_IO_COVER_WRITER_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Write a clustering to a file.
 */

class CoverWriter final {
public:
    void write(Cover &zeta, const std::string &path) const;
};
} // namespace NetworKit

#endif // NETWORKIT_IO_COVER_WRITER_HPP_
