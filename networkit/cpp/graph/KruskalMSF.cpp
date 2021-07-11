// no-networkit-format
/*
 * KruskalMSF.cpp
 *
 *  Created on: 03.09.2015
 *      Author: Henning
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/KruskalMSF.hpp>
#include <networkit/graph/SpanningForest.hpp>
#include <networkit/structures/UnionFind.hpp>

namespace NetworKit {

struct MyEdge {
    node from;
    node to;
    edgeweight weight;

    MyEdge(node u, node v, edgeweight ew) {
        from = u;
        to = v;
        weight = ew;
    }

    MyEdge() {
        from = none;
        to = none;
        weight = 0;
    }

    bool operator<(const MyEdge other) const {
        return this->weight > other.weight; // Note the switch in the operator!
    }
};


KruskalMSF::KruskalMSF(const Graph& G): SpanningForest(G) {}

void KruskalMSF::run() {
    if (true || G->isWeighted()) { // FIXME: remove true when SpanningForest is fixed!
        count z = G->upperNodeIdBound();
        forest = GraphTools::copyNodes(*G);
        UnionFind uf(z);

        // sort edges in decreasing weight order
        std::vector<MyEdge> sortedEdges; // (m);
        G->forEdges([&](node u, node v, edgeweight ew, edgeid) {
            MyEdge myEdge(u, v, ew);
            sortedEdges.push_back(myEdge);
        });
        Aux::Parallel::sort(sortedEdges.begin(), sortedEdges.end());

        // process in decreasing weight order
        for (auto e: sortedEdges) {
            node u = e.from;
            node v = e.to;
            INFO("process edge (", u, ", ", v, ") with weight ", e.weight);
            assert(u < z);
            assert(v < z);

            // if edge does not close cycle, add it to tree
            if (uf.find(u) != uf.find(v)) {
                forest.addEdge(u, v);
                uf.merge(u, v);
            }
        }
    }
    else {
        SpanningForest sf(*G);
        sf.run();
        forest = sf.getForest();
    }
}

} /* namespace NetworKit */
