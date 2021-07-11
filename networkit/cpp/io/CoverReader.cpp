
#include <networkit/io/CoverReader.hpp>

#include <fstream>
#include <sstream>

namespace NetworKit {

Cover CoverReader::read(const std::string &path, Graph &G) {
    std::ifstream file;
    file.open(path);
    if (!file.good()) {
        throw std::runtime_error("unable to read from file");
    }
    Cover communities(G.upperNodeIdBound());
    std::string line;
    count i = 0;
    node current;

    while (std::getline(file, line)) {
        if (line.substr(0, 1) != "#") {
            communities.setUpperBound(i + 1);
            std::stringstream linestream(line);
            while (linestream >> current) {
                communities.addToSubset(i, current);
            }
            ++i;
        }
    }

    file.close();
    return communities;
}

} // namespace NetworKit
