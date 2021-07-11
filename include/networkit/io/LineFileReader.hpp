/*
 * LineFileReader.hpp
 *
 *  Created on: 04.10.2013
 *      Author: cls
 */

#ifndef NETWORKIT_IO_LINE_FILE_READER_HPP_
#define NETWORKIT_IO_LINE_FILE_READER_HPP_

#include <fstream>
#include <string>
#include <vector>

namespace NetworKit {

/**
 * @ingroup io
 * Reads a file and puts each line as a string into a vector
 */
class LineFileReader final {
public:
    std::vector<std::string> read(std::string path) {
        std::ifstream file;
        std::string line; // the current line
        file.open(path);

        std::vector<std::string> data;

        while (file.good()) {
            std::getline(file, line);
            data.push_back(line);
        }

        file.close();

        return data;
    }
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_LINE_FILE_READER_HPP_
