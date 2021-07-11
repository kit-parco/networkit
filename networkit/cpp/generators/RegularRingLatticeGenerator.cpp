// no-networkit-format
/*
* RegularRingLatticeGenerator.cpp
*
*  Created on: 09.07.2014
*      Author: Simon Bischof
*/

#include <networkit/generators/RegularRingLatticeGenerator.hpp>

namespace NetworKit {

RegularRingLatticeGenerator::RegularRingLatticeGenerator(count nNodes, count numberOfNeighbors)
    : nNodes(nNodes), nNeighbors(std::min(numberOfNeighbors, nNodes / 2 - 1)) {}

Graph RegularRingLatticeGenerator::generate() {
    Graph G(nNodes);
    for (count i = 0; i < nNodes; i++) {
        for (count j = 1; j <= nNeighbors; j++) {
            G.addEdge(i, (i + j) % nNodes);
        }
    }

    return G;
}

} /* namespace NetworKit */
