/*
 * METISGraphReader.hpp
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_METIS_GRAPH_READER_HPP_
#define NETWORKIT_IO_METIS_GRAPH_READER_HPP_

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Reader for the METIS file format documented in [1]
 *
 * [1] http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf
 */
class METISGraphReader final : public GraphReader {
public:
    METISGraphReader() = default;

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
#endif // NETWORKIT_IO_METIS_GRAPH_READER_HPP_
