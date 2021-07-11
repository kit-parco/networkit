// no-networkit-format
/*
 * PartitionReader.cpp
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt
 */

#include <networkit/io/PartitionReader.hpp>

namespace NetworKit {

Partition PartitionReader::read(const std::string &path) {

    std::ifstream file(path);

    // check if file readable
    if (!file) {
        throw std::runtime_error("invalid clustering file");
    }


    Partition zeta(0);

    std::string line;
    index omega = 0;
    while(std::getline(file, line)) {
        if (line.substr(0, 1) == "*" || line.substr(0, 1) == "#") continue;

        index c = std::atoi(line.c_str());
        // extend the partition by one entry and store the cluster id
        zeta[zeta.extend()] = c;
        if (c != none) {
            omega = std::max(c, omega);
        }
    }

    zeta.setUpperBound(omega+1);

    return zeta;
}

} /* namespace NetworKit */
