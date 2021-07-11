// no-networkit-format
/*
 * MatchingCoarsening.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt
 */

#include <networkit/coarsening/MatchingCoarsening.hpp>

namespace NetworKit {

MatchingCoarsening::MatchingCoarsening(const Graph& G, const Matching& M, bool noSelfLoops) : GraphCoarsening(G), M(M), noSelfLoops(noSelfLoops) {
    if (G.isDirected()) throw std::runtime_error("Only defined for undirected graphs.");
}

void MatchingCoarsening::run() {
    count n = G->numberOfNodes();
    index z = G->upperNodeIdBound();
    count cn = n - M.size(*G);
    Graph cG(cn, true);

    // compute map: old ID -> new coarse ID
    index idx = 0;
    std::vector<node> mapFineToCoarse(z, none);
    G->forNodes([&](node v) { // TODO: difficult in parallel
        index mate = M.mate(v);
        if (mate == v) DEBUG("Node ", v, " is its own matching!");
        assert(mate != v);
        if ((mate == none) || (v < mate)) {
            // vertex v is carried over to the new level
            mapFineToCoarse[v] = idx;
            ++idx;
        }
        else {
            // vertex v is not carried over, receives ID of mate
            mapFineToCoarse[v] = mapFineToCoarse[mate];
        }
        assert(mapFineToCoarse[v] != none);
        assert(mapFineToCoarse[v] < cn);
    });

    G->forNodes([&](node v) { // TODO: difficult in parallel
        G->forNeighborsOf(v, [&](node u, edgeweight ew) {
            node cv = mapFineToCoarse[v];
            node cu = mapFineToCoarse[u];
            if ((v <= u) && (! noSelfLoops || (cv != cu))) {
                cG.increaseWeight(cv, cu, ew);
            }
        });
    });

    Gcoarsened = std::move(cG);
    nodeMapping = std::move(mapFineToCoarse);

    hasRun = true;
}

} /* namespace NetworKit */
