// no-networkit-format
/*
* DorogovtsevMendesGenerator.cpp
*
*  Created on: 27.05.2014
*      Author: Christian Staudt
*/

#include <tuple>

#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

DorogovtsevMendesGenerator::DorogovtsevMendesGenerator(count nNodes): nNodes(nNodes) {}

Graph DorogovtsevMendesGenerator::generate() {
    Graph G;
    // create initial triangle
    node s1 = G.addNode();
    node s2 = G.addNode();
    node s3 = G.addNode();
    G.addEdge(s1, s2);
    G.addEdge(s2, s3);
    G.addEdge(s3, s1);

    for (count i = 0; i < (nNodes - 3); ++i) {
        node u;
        node v;
        std::tie(u, v) = GraphTools::randomEdge(G);
        node w = G.addNode();
        G.addEdge(w, u);
        G.addEdge(w, v);
    }

    G.shrinkToFit();
    return G;
}

} /* namespace NetworKit */
