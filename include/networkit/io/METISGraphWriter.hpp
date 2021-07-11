/*
 * METISGraphWriter.hpp
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_METIS_GRAPH_WRITER_HPP_
#define NETWORKIT_IO_METIS_GRAPH_WRITER_HPP_

#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class METISGraphWriter final : public GraphWriter {

public:
    METISGraphWriter() = default;

    void write(const Graph &G, const std::string &path) override;

    void write(const Graph &G, bool weighted, const std::string &path);
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_METIS_GRAPH_WRITER_HPP_
