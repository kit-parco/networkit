/*
 * METISParser.hpp
 *
 *  Created on: 27.11.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_METIS_PARSER_HPP_
#define NETWORKIT_IO_METIS_PARSER_HPP_

#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <networkit/Globals.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Parser for the METIS file format.
 */
class METISParser final {

    std::ifstream graphFile;

public:
    /**
     * Constructor for the METIS Parser.
     *
     * @param[in]  path  file path
     */
    METISParser(const std::string &path);

    /**
     * Get the METIS graph file header
     */
    std::tuple<count, count, index, count> getHeader();

    /**
     * Test if graph file has a next line.
     */
    bool hasNext();

    /**
     * Get adjacencies from the next line in the METIS graph file.
     *
     * @param[in]  ignoreFirst  number of values to ignore [in case the METIS file contains node
     * weightes]
     */
    std::vector<node> getNext(count ignoreFirst = 0);

    /**
     * Get adjacencies with edge weights from the next line in the METIS graph file.
     *
     * @param[in]  ignoreFirst  number of values to ignore [in case the METIS file contains node
     * weightes]
     */
    std::vector<std::pair<node, double>> getNextWithWeights(count ignoreFirst = 0);
};
} /* namespace NetworKit */
#endif // NETWORKIT_IO_METIS_PARSER_HPP_
