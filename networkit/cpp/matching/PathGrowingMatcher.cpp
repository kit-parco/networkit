// no-networkit-format
/*
 * PathGrowingMatcher.cpp
 *
 *  Created on: Jun 13, 2013
 *      Author: Henning
 */

#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/matching/PathGrowingMatcher.hpp>

namespace NetworKit {

PathGrowingMatcher::PathGrowingMatcher(const Graph& G): Matcher(G) { checkInput(); }

PathGrowingMatcher::PathGrowingMatcher(const Graph &G, const std::vector<double> &edgeScores)
    : Matcher(G, edgeScores) {
    checkInput();
}

void PathGrowingMatcher::checkInput() const {
    if (G->isDirected())
        throw std::invalid_argument(
            "The input graph is directed, this algorithm only supports undirected graphs.");

    if (G->numberOfSelfLoops() > 0)
        throw std::invalid_argument("The input graph has self-loops, this algorithm only supports "
                                    "graphs without self-loops. ");
}

void PathGrowingMatcher::run() {
    count z = G->upperNodeIdBound();

    // init matching to empty
    Matching m1(z);
    Matching m2(z);
    bool takeM1 = true;

    // PQ to retrieve vertices with degree > 0 quickly
    int64_t minKey = -((int64_t) G->numberOfNodes());
    int64_t maxKey = 0;
    Aux::BucketPQ bpq(z, minKey, maxKey);

    // degrees tracks degree of vertices,
    // avoids to make a copy of the graph and
    // delete vertices and edges explicitly.
    std::vector<int64_t> degrees(z, 0);

    // alive tracks if vertices are alive or not in the algorithm
    std::vector<bool> alive(z, false);

    G->forNodes([&](node u) {
        degrees[u] = G->degree(u);
        alive[u] = (degrees[u] > 0);
        if (alive[u]) {
            bpq.insert(-degrees[u], u); // minus to get extractMin functionality
        }
    });

    count numEdges = G->numberOfEdges();

    // main loop
    while (numEdges > 0) {
        // use vertex with positive degree
        std::pair<int64_t, node> mini = bpq.extractMin();
        node v = mini.second;
        assert(v != none);

        // path growing
        while (degrees[v] > 0) {
            // find heaviest incident edge
            node bestNeighbor = 0;
            edgeweight bestWeight = 0;

            if (edgeScoresAsWeights) {
                G->forEdgesOf(v, [&](node, node u, edgeid eid) {
                    if (alive.at(u)) {
                        if (edgeScores.at(eid) > bestWeight) {
                            bestNeighbor = u;
                            bestWeight = edgeScores.at(eid);
                        }
                    }
                });
            } else {
                G->forEdgesOf(v, [&](node, node u, edgeweight weight) {
                    if (alive.at(u)) {
                        if (weight > bestWeight) {
                            bestNeighbor = u;
                            bestWeight = weight;
                        }
                    }
                });
            }


            if (takeM1) {
                // add edge to m1
                m1.match(v, bestNeighbor);
                takeM1 = false;
            }
            else {
                // add edge to m2
                m2.match(v, bestNeighbor);
                takeM1 = true;
            }

            // remove current vertex and its incident edges from graph
            G->forEdgesOf(v, [&](node, node u) {
                if (alive.at(u)) {
                    degrees.at(u)--;
                    numEdges--;

                    if (degrees.at(u) == 0) {
                        // singleton node can be removed
                        bpq.remove(u);
                        alive[u] = false;
                    }
                    else {
                        bpq.changeKey(-degrees.at(u), u);
                    }
                }
            });
            alive[v] = false;
            bpq.remove(v);

            // start next iteration from best neighbor
            v = bestNeighbor;
        }
    }

    // return the heavier one of the two
    edgeweight weight1 {0};
    if (edgeScoresAsWeights) {
        G->forEdges([&](node, node, edgeid eid){
            weight1 += edgeScores.at(eid);
        });
    } else {
        weight1 = m1.weight(*G);
    }
    edgeweight weight2 = 0.;
    if (edgeScoresAsWeights) {
        G->forEdges([&](node, node, edgeid eid){
            weight2 += edgeScores.at(eid);
        });
    } else {
        weight2 = m2.weight(*G);
    }
    INFO("weight of first matching: ", weight1);
    INFO("weight of second matching: ", weight2);

    if (weight1 > weight2)
        M = m1;
    else
        M = m2;

    hasRun = true;
}

} /* namespace NetworKit */
