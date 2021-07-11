
#ifndef NETWORKIT_IO_DOT_PARTITION_WRITER_HPP_
#define NETWORKIT_IO_DOT_PARTITION_WRITER_HPP_

#include <map>

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class DotPartitionWriter final {
public:
    void write(Graph &graph, Partition &zeta, const std::string &path) const;

    std::map<index, double> createHueMap(Graph &graph, Partition &zeta) const;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_DOT_PARTITION_WRITER_HPP_
