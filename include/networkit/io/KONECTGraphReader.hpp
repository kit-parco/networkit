/*
 * KONECTGraphReader.hpp
 *
 *  Created on: 11.05.2018
 *      Author: Roman Bange
 */

#ifndef NETWORKIT_IO_KONECT_GRAPH_READER_HPP_
#define NETWORKIT_IO_KONECT_GRAPH_READER_HPP_

#include <unordered_map>

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphReader.hpp>

namespace NetworKit {
class KONECTGraphReader final : public GraphReader {

public:
    /*
     * If the input graph has multiple edges, you can specify on how these edges are handled.
     * Keep in mind that NetworKit node id's start with 0 while most KONECT graphs start with 1.
     * See GraphReader.h for a closer description of the parameters.
     *
     * @param[in]  remapNodes  specifies whether node ids should be remapped if non consecutive
     * @param[in]  handlingmethod  specifies how multiple edges should be handled (only relevant if
     * graph with multiple edges is given)
     */
    KONECTGraphReader(bool remapNodes = false,
                      MultipleEdgesHandling handlingmethod = DISCARD_EDGES);

    /**
     * Given the path of an input file, read the graph contained.
     *
     * @param[in]  path  input file path
     */
    Graph read(const std::string &path) override;

private:
    bool remapNodes;
    MultipleEdgesHandling multipleEdgesHandlingMethod;
};
} /* namespace NetworKit */
#endif // NETWORKIT_IO_KONECT_GRAPH_READER_HPP_
