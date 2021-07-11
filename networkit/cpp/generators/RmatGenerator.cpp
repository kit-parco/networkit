// no-networkit-format
/*
 * RmatGenerator.cpp
 *
 *  Created on: 18.03.2014
 *      Author: Henning, cls
 */

#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/generators/RmatGenerator.hpp>

namespace NetworKit {

RmatGenerator::RmatGenerator(count scale, count edgeFactor, double a, double b, double c, double d, bool weighted, count reduceNodes):
    scale(scale), edgeFactor(edgeFactor), a(a), b(b), c(c), weighted(weighted), reduceNodes(reduceNodes)
{
    if (scale > 63) throw std::runtime_error("Cannot generate more than 2^63 nodes");
    double sum = a+b+c+d;
    INFO("sum of probabilities: ", sum);
    if (!Aux::NumericTools::equal(sum, 1.0, 0.0001)) throw std::runtime_error("Probabilities in Rmat have to sum to 1.");
    defaultEdgeWeight = 1.0;
}

Graph RmatGenerator::generate() {
    double n = (1 << scale);
    if (n <= reduceNodes) {
        throw std::runtime_error("Error, shall delete more nodes than the graph originally has");
    }
    // when nodes are deleted, all nodes have less neighbors
    count numEdges = n * edgeFactor * n * 1.0 / static_cast<double>(n - reduceNodes);
    count wantedEdges = (n - reduceNodes) * edgeFactor;
    Graph G(n - reduceNodes, weighted);
    double ab = a+b;
    double abc = ab+c;

    auto quadrant([&]() {
        double r = Aux::Random::probability();
        TRACE("r: ", r);

        if (r <= a) {
            return 0;
        }
        else if (r <= ab) {
            return 1;
        }
        else if (r <= abc) {
            return 2;
        }
        else return 3;
    });

    auto drawEdge([&]() {
        node u = 0;
        node v = 0;
        for (index i = 0; i < scale; ++i) {
            count q = quadrant();
            u = u << 1;
            v = v << 1;
            u = u | (q >> 1);
            v = v | (q & 1);
        }

        return std::make_pair(u, v);
    });

    if (reduceNodes > 0) {
        std::vector<node> nodemap(n, 0);

        for (count deletedNodes = 0; deletedNodes < reduceNodes;) {
            node u = Aux::Random::index(n);
            if (nodemap[u] == 0) {
                nodemap[u] = none;
                ++deletedNodes;
            }
        }

        for (node i = 0, u = 0; i < n; ++i) {
            if (nodemap[i] == 0) {
                nodemap[i] = u;
                ++u;
            }
        }

        node u, v;
        if (weighted) {
            for (index e = 0; e < numEdges; ++e) {
                std::tie(u, v) = drawEdge();
                u = nodemap[u];
                v = nodemap[v];
                if (u != none && v != none) {
                    G.increaseWeight(u, v, defaultEdgeWeight);
                }
            }
        } else {
            while (G.numberOfEdges() < wantedEdges) {
                std::tie(u, v) = drawEdge();
                u = nodemap[u];
                v = nodemap[v];
                if (u != none && v != none && u != v && !G.hasEdge(u, v)) {
                    G.addEdge(u, v);
                }
            }
        }
    } else {
        if (weighted) {
            for (index e = 0; e < numEdges; ++e) {
                std::pair<node, node> drawnEdge = drawEdge();
                G.increaseWeight(drawnEdge.first, drawnEdge.second, defaultEdgeWeight);
            }
        } else {
            while (G.numberOfEdges() < wantedEdges) {
                std::pair<node, node> drawnEdge = drawEdge();
                if (!G.hasEdge(drawnEdge.first, drawnEdge.second)) {
                    G.addEdge(drawnEdge.first, drawnEdge.second);
                }
            }
        }
    }

    G.shrinkToFit();
    return G;
}

} /* namespace NetworKit */
