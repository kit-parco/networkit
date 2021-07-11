/*
 * WattsStrogatzGenerator.cpp
 *
 *  Created on: 09.07.2014
 *      Author: Simon Bischof
 */

#include <limits>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/RegularRingLatticeGenerator.hpp>
#include <networkit/generators/WattsStrogatzGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

WattsStrogatzGenerator::WattsStrogatzGenerator(count nNodes, count numberOfNeighbors, double p)
    : nNodes(nNodes), nNeighbors(numberOfNeighbors), p(p) {
    if (numberOfNeighbors * 2 >= nNodes - 1) {
        throw std::runtime_error("nNeighbors*2 cannot be equal to nNodes-1.");
    }
    if (numberOfNeighbors > nNodes / 2 - 1) {
        nNeighbors = nNodes / 2 - 1;
    }
}

Graph WattsStrogatzGenerator::generate() {
    // generate regular ring lattice as initial graph
    auto G = RegularRingLatticeGenerator(nNodes, nNeighbors).generate();
    std::uniform_int_distribution<node> dist(0, nNodes - 1);
    std::mt19937_64 generator(Aux::Random::integer());
    std::uniform_real_distribution<double> prob{0.0, std::nexttoward(1.0, 2.0)};

    // list all valid end nodes for rewiring an edge
    auto validEndNode = [&G](const node u, const node v) { return u != v && !G.hasEdge(u, v); };

    std::vector<node> prevNeighbors;
    prevNeighbors.resize(GraphTools::maxDegree(G));

    // rewire according to Watts and Strogatz model
    for (const auto u : G.nodeRange()) {
        if (G.degree(u) == (nNodes - 1))
            continue;
        /* save edges before rewiring incident edges for a node
         * because they may get changed through rewiring
         */
        count prevNeighborsSize = 0;
        G.forNeighborsOf(u, [&](const node v) {
            if (u < v) {
                prevNeighbors[prevNeighborsSize++] = v;
            }
        });
        std::for_each(prevNeighbors.begin(), prevNeighbors.begin() + prevNeighborsSize,
                      [&](const node v) {
                          if (prob(generator) < p) {
                              node w = none;
                              do {
                                  w = dist(generator);
                              } while (!validEndNode(u, w));
                              G.removeEdge(u, v);
                              G.addEdge(u, w);
                          }
                      });
    }

    return G;
}

} /* namespace NetworKit */
