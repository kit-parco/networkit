// no-networkit-format
/*
 *
 */

#include <stdexcept>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/centrality/PermanenceCentrality.hpp>

NetworKit::PermanenceCentrality::PermanenceCentrality(const NetworKit::Graph &G, const NetworKit::Partition &P): Algorithm(), G(G), P(P) {
}

void NetworKit::PermanenceCentrality::run() {
    Aux::SignalHandler handler;

    auto isOutEdge = [&](node u, node v) {
        return G.degree(u) > G.degree(v) || (G.degree(u) == G.degree(v) && u < v);
    };

    inEdges.resize(G.numberOfEdges());
    {
        // bucket sort of nodes by degree
        count n = G.numberOfNodes();
        std::vector<node> sortedNodes(n);
        {
            std::vector<index> nodePos(n, 0);

            G.forNodes([&](node u) {
                ++nodePos[G.degree(u)];
            });

            // exclusive prefix sum
            index tmp = nodePos[0];
            index sum = tmp;
            nodePos[0] = 0;

            for (index i = 1; i < nodePos.size(); ++i) {
                tmp = nodePos[i];
                nodePos[i] = sum;
                sum += tmp;
            }

            G.forNodes([&](node u) {
                sortedNodes[nodePos[G.degree(u)]++] = u;
            });
        }

        INFO("Sorted nodes");

        handler.assureRunning();

        // count in-degrees
        std::vector<count> inDeg(G.upperNodeIdBound() + 1);

        G.forEdges([&](node u, node v) {
            if (isOutEdge(u, v)) {
                inDeg[v]++;
            } else {
                inDeg[u]++;
            }
        });

        // calculate positions of incoming edges in the adjacency array
        inBegin.swap(inDeg);

        {
            // exclusive prefix sum
            index tmp = inBegin[0];
            index sum = tmp;
            inBegin[0] = 0;

            for (index i = 1; i < inBegin.size(); ++i) {
                tmp = inBegin[i];
                inBegin[i] = sum;
                sum += tmp;
            }

            inBegin[G.upperNodeIdBound()] = sum;
        }

        handler.assureRunning();

        INFO("Got edge positions");

        // the huge shuffle that moves all edges into their position
        for (node u : sortedNodes) {
            G.forEdgesOf(u, [&](node v) {
                if (isOutEdge(u, v)) {
                    inEdges[inBegin[v]++] = u;
                }
            });
        }

        INFO("Created edge array");

        // begin became end, so move it to be begin again!
        index tmp = 0;
        for (index i = 0; i < inBegin.size(); ++i) {
            std::swap(inBegin[i], tmp);
        }
    }

    marker.clear();
    marker.resize(G.upperNodeIdBound(), false);

    INFO("Finished copying G");

    hasRun = true;
}

double NetworKit::PermanenceCentrality::getIntraClustering(NetworKit::node u) {
    assureFinished();
    count numNeighbors = 0;
    index C = P[u];

    G.forNeighborsOf(u, [&](node y) {
        marker[y] = (P[y] == C);
        numNeighbors += marker[y];
    });

    count numTriangles = 0;

    G.forNeighborsOf(u, [&](node y) {
        if (marker[y]) {
            for (index i = inBegin[y]; i < inBegin[y + 1]; ++i) {
                node z = inEdges[i];
                numTriangles += marker[z];
            }
        }
    });

    G.forNeighborsOf(u, [&](node y) {
        marker[y] = false;
    });

    if (numNeighbors < 2) return 0;

    return numTriangles * 1.0 / (0.5 * static_cast<double>(numNeighbors * (numNeighbors - 1))); // triangles are counted only once, so divide by 2 here!
}

double NetworKit::PermanenceCentrality::getPermanence(NetworKit::node u) {
    std::map<index, count> strength;
    index C = P[u];
    G.forNeighborsOf(u, [&](node y) {
        ++strength[P[y]];
    });

    count maxNeighborStrength = 0;
    for (auto n_s : strength) {
        if (n_s.first != C && n_s.second > maxNeighborStrength) {
            maxNeighborStrength = n_s.second;
        }
    }

    //  http://dl.acm.org/citation.cfm?doid=2623330.2623707: If the vertex has no external connections, F1 is just the value of the inter- nal connections.
    if (maxNeighborStrength == 0) maxNeighborStrength = 1;

    return strength[C] * 1.0 / maxNeighborStrength * 1.0 /G.degree(u) - (1.0 - getIntraClustering(u));
}
