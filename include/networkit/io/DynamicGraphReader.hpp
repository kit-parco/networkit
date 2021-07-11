/*
 * DynamicGraphReader.hpp
 *
 *  Created on: 01.06.2013
 *      Author: cls
 */

#ifndef NETWORKIT_IO_DYNAMIC_GRAPH_READER_HPP_
#define NETWORKIT_IO_DYNAMIC_GRAPH_READER_HPP_

#include <networkit/dynamics/GraphEventProxy.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class DynamicGraphReader {

public:
    /**
     * @param[in] path path to dynamic graph file
     * @param[in] Gproxy graph event proxy receives the events from the file
     */
    virtual void read(std::string path, GraphEventProxy &Gproxy) = 0;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_DYNAMIC_GRAPH_READER_HPP_
