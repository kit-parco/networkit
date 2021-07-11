// no-networkit-format
/*
 * EdgeListPartitionReader.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: forigem
 */

#include <networkit/io/EdgeListPartitionReader.hpp>

namespace NetworKit {

EdgeListPartitionReader::EdgeListPartitionReader(node firstNode, char sepChar) : firstNode(firstNode), sepChar(sepChar) {
}

Partition EdgeListPartitionReader::read(const std::string &path) {
    std::ifstream file(path);

    // check if file readable
    if (!file) {
        throw std::runtime_error("invalid clustering file");
    }

    std::string line;
    index newOmega = 0;
    Partition zeta(0); // start without nodes
    while(std::getline(file, line)) {
        std::vector<std::string> split = Aux::StringTools::split(line, sepChar);
        if (split[0] != "#") {
            index c = std::atoi(split[1].c_str());
            // NetworKit uses zero-based node ids, adapt input accordingly
            index v = std::atoi(split[0].c_str()) - firstNode;
            // add nodes until we reach the current node
            while (v >= zeta.numberOfElements()) {
                zeta.extend();
            }
            // calculate the upper bound of the cluster ids
            newOmega = std::max(newOmega, c);
            zeta[v] = c;
        }
    }

    zeta.setUpperBound(newOmega+1);
    return zeta;
}

} /* namespace NetworKit */
