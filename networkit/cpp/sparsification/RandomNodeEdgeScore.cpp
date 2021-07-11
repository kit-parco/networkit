// no-networkit-format
/*
 * RandomNodeEdgeScore.cpp
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/sparsification/RandomNodeEdgeScore.hpp>

namespace NetworKit {

RandomNodeEdgeScore::RandomNodeEdgeScore(const Graph& G, double rneRatio) : EdgeScore<double>(G), rneRatio(rneRatio) {
}

void RandomNodeEdgeScore::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    Graph sparseGraph = *G;
    std::vector<double> workScores(G->upperEdgeIdBound(), 0);
    count numRemoved = 0;
    std::vector< std::pair<node, node> > uniformlyRandomEdges;

    while (sparseGraph.numberOfEdges() > 0) {
        if (Aux::Random::real() >= rneRatio) { // uniformly random
            bool edgeFound = false;

            while (!edgeFound) {
                if (uniformlyRandomEdges.empty()) {
                    uniformlyRandomEdges = GraphTools::randomEdges(sparseGraph, sparseGraph.numberOfEdges() * (1.0 - rneRatio) + 20);
                }

                auto edge = uniformlyRandomEdges.back();
                uniformlyRandomEdges.pop_back();

                if (sparseGraph.hasEdge(edge.first, edge.second)) {
                    edgeid id = sparseGraph.edgeId(edge.first, edge.second);

                    workScores[id] = numRemoved * 1.0 / G->numberOfEdges();

                    sparseGraph.removeEdge(edge.first, edge.second);

                    edgeFound = true;
                    ++numRemoved;
                }
            }
        } else { // random node - edge
            const auto edge = GraphTools::randomEdge(sparseGraph);

            edgeid id = sparseGraph.edgeId(edge.first, edge.second);

            workScores[id] = numRemoved * 1.0 / G->numberOfEdges();

            sparseGraph.removeEdge(edge.first, edge.second);

            ++numRemoved;
        }
    }

    scoreData = std::move(workScores);
    hasRun = true;
}

double RandomNodeEdgeScore::score(node, node) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

double RandomNodeEdgeScore::score(edgeid) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
