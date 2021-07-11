// no-networkit-format
/*
 * EdgeCut.cpp
 *
 *  Created on: Jun 20, 2013
 *      Author: Henning
 */

#include <networkit/community/EdgeCut.hpp>

namespace NetworKit {

double EdgeCut::getQuality(const Partition& zeta, const Graph& G) {
    double cutWeight = 0.0;
    G.forEdges([&](node u, node v, edgeweight w) {
        if (zeta[u] != zeta[v]) {
            cutWeight += w;
        }
    });
    return cutWeight;
}

} /* namespace NetworKit */
