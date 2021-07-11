/*
 * GraphIO.cpp
 *
 *  Created on: 09.01.2013
 *      Author: Christian Staudt
 */

#include <fstream>
#include <networkit/io/GraphIO.hpp>

namespace NetworKit {

void GraphIO::writeEdgeList(const Graph &G, const std::string &path) {

    std::ofstream file;
    file.open(path.c_str());

    G.forEdges([&](const node v, const node w) { file << v << " " << w << '\n'; });

    file.close();
}

void GraphIO::writeAdjacencyList(const Graph &G, const std::string &path) {
    std::ofstream file;
    file.open(path.c_str());

    G.forNodes([&](const node v) {
        file << v;
        G.forNeighborsOf(v, [&](const node x) { file << " " << x; });
        file << '\n';
    });
}

} /* namespace NetworKit */
