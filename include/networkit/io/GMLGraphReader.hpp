/*
 * GMLGraphReader.hpp
 *
 *  Created on: 18.09.2014
 *      Author: Maximilian Vogel (maximilian.vogel@student.kit.edu)
 */

#ifndef NETWORKIT_IO_GML_GRAPH_READER_HPP_
#define NETWORKIT_IO_GML_GRAPH_READER_HPP_

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class GMLGraphReader final : public GraphReader {
public:
    GMLGraphReader() = default;

    /**
     * Takes a file path as parameter and returns a graph file.
     *
     * @param[in]  path  file path
     *
     * @param[out]  the graph read from file
     */
    Graph read(const std::string &path) override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_GML_GRAPH_READER_HPP_
