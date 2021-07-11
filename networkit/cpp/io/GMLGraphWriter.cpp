/*
 * GMLGraphWriter.cpp
 *
 *  Created on: 05.08.2013
 *      Author: Stefan Bertsch
 */

#include <fstream>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/io/GMLGraphWriter.hpp>

namespace NetworKit {

void GMLGraphWriter::write(const Graph &G, const std::string &path) {
    std::ofstream file(path);
    Aux::enforceOpened(file);

    file << "graph [\n";
    if (G.isDirected()) {
        file << "  directed 1\n";
    }

    G.forNodes([&](node u) {
        file << "  node [\n";
        file << "    id " << u << "\n";
        file << "  ]\n";
    });

    G.forEdges([&](node u, node v) {
        file << "  edge [\n";
        file << "    source " << u << "\n";
        file << "    target " << v << "\n";
        file << "  ]\n";
    });
    file << "]\n";
}

} /* namespace NetworKit */
