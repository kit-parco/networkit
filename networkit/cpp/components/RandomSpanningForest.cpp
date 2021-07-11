// no-networkit-format
/*
 * RandomSpanningForest.cpp
 *
 *  Created on: 06.09.2015
 *      Author: Henning
 */

#include <unordered_set>

#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/components/RandomSpanningForest.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

RandomSpanningForest::RandomSpanningForest(const Graph& G):
        SpanningForest(G) {}

void RandomSpanningForest::run() {
    // handle disconnected graphs:
    // determine connected components first
    // then start random walk in each component!
    ConnectedComponents cc(*G);
    cc.run();
    const auto comps = cc.getComponents();

    forest = GraphTools::copyNodes(*G);
    for (const auto &comp: comps) {
        std::unordered_set<node> visited;

        // find and process random root
        index rand = Aux::Random::integer(comp.size() - 1);
        node curr = comp[rand];
        visited.insert(curr);

        // random walk starting from root
        while (visited.size() < comp.size()) {
            // get random neighbor
            node neigh = GraphTools::randomNeighbor(*G, curr);

            // if not seen before, insert tree edge
            if (visited.count(neigh) == 0) {
                forest.addEdge(curr, neigh);
                visited.insert(neigh);
            }

            // move to neighbor
            curr = neigh;
        }
    }
}

} /* namespace NetworKit */
