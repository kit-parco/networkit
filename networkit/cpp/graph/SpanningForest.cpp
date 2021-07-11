// no-networkit-format
/*
 * SpanningForest.cpp
 *
 *  Created on: 06.09.2015
 *      Author: Henning
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/SpanningForest.hpp>

namespace NetworKit {

void SpanningForest::run() {
    forest = GraphTools::copyNodes(*G);
    std::vector<bool> visited(G->upperNodeIdBound(), false);

    G->forNodes([&](node s){
        if (! visited[s]) {
            Traversal::BFSEdgesFrom(*G, s, [&](node u, node v, edgeweight w, edgeid) {
                visited[u] = true;
                visited[v] = true;
                forest.addEdge(u, v, w);
            });
        }
    });

    INFO("tree edges in SpanningForest: ", forest.numberOfEdges());
}

} /* namespace NetworKit */
