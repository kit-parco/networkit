/*
 * EdgeListReader.hpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef NETWORKIT_IO_EDGE_LIST_READER_HPP_
#define NETWORKIT_IO_EDGE_LIST_READER_HPP_

#include <map>
#include <string>

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * A reader for various edge list formats, in which each line contains an edge as
 * two node ids.
 *
 */
class EdgeListReader final : public GraphReader {

public:
    EdgeListReader() = default; // nullary constructor for Python shell

    /**
     * @param[in]  separator  character used to separate nodes in an edge line
     * @param[in]  firstNode  index of the first node in the file
     * @param[in]  commentPrefix  prefix of comment lines
     * @param[in]  continuous  boolean to specify if node ids are continuous
     * @param[in]  directed  read graph as directed
     */
    EdgeListReader(char separator, node firstNode, const std::string &commentPrefix = "#",
                   bool continuous = true, bool directed = false);

    /**
     * Given the path of an input file, read the graph contained.
     *
     * @param[in]  path  input file path
     */
    Graph read(const std::string &path) override;

    /**
     * Return the node map, in case node ids are not continuous
     */
    const std::map<std::string, node> &getNodeMap() const;

private:
    char separator; //!< character separating nodes in an edge line
    std::string commentPrefix;
    node firstNode;
    bool continuous;
    std::map<std::string, node> mapNodeIds;
    bool directed;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_EDGE_LIST_READER_HPP_
