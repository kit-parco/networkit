/*
 * SNAPGraphReader.hpp
 *
 *  Created on: 04.05.2018
 *      Author: Alexander van der Grinten
 */

#ifndef NETWORKIT_IO_SNAP_GRAPH_READER_HPP_
#define NETWORKIT_IO_SNAP_GRAPH_READER_HPP_

#include <unordered_map>

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class SNAPGraphReader final : public GraphReader {
    std::unordered_map<node, node> nodeIdMap;
    bool directed;
    count nodeCount;
    bool remapNodes;

public:
    /**
     * Default constructor for the SNAPGraphReader.
     * NOTE: Keep in mind that many SNAP graphs do not have consecutive node ids.
     * Passing a "false" in remapNodes forces the reader to add nodes for the highest occuring node
     * id. This may result in a lot of isolated nodes, but further increases the speed of x2(in
     * benchmarking)
     *
     * @param[in]  directed  reads in the graph as directed, if set to true
     * @param[in]  remapNodes  indicates whether nodes should be remapped to other node ids in order
     * to create consecutive node ids
     * @param[in]  nodeCount  nodeCount is used to preallocated memory for the number of nodes
     */
    SNAPGraphReader(bool directed = false, const bool &remapNodes = true,
                    const count &nodeCount = 0);

    /**
     * Given the path of an input file, read the graph contained.
     *
     * @param[in]  path  input file path
     */
    Graph read(const std::string &path) override;
};

} /* namespace NetworKit */

#endif // NETWORKIT_IO_SNAP_GRAPH_READER_HPP_
