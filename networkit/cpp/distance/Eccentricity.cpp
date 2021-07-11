// no-networkit-format
/*
 * Eccentricity.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include <networkit/distance/Eccentricity.hpp>
#include <networkit/graph/BFS.hpp>

namespace NetworKit {

std::pair<node, count> Eccentricity::getValue(const Graph& G, node u) {
    assert(G.hasNode(u));
    count ecc = 0;
    node res = none;
    Traversal::BFSfrom(G, u, [&](node v, count dist) {
        ecc = dist;
        res = v;
    });
    assert(res != none);
    return {res, ecc}; // pair.first is argmax node
}

} /* namespace NetworKit */
