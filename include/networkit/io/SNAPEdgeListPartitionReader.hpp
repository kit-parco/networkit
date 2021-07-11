/*
 * SNAPEdgeListPartitionReader.hpp
 *
 *  Created on: Jun 20, 2013
 *      Author: forigem
 */

#ifndef NETWORKIT_IO_SNAP_EDGE_LIST_PARTITION_READER_HPP_
#define NETWORKIT_IO_SNAP_EDGE_LIST_PARTITION_READER_HPP_

#include <unordered_map>

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Reads the clustering files from the SNAP collection.
 */
class SNAPEdgeListPartitionReader final {
public:
    Cover read(const std::string &path, std::unordered_map<node, node> &mapNodeIds, Graph &G);
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_SNAP_EDGE_LIST_PARTITION_READER_HPP_
