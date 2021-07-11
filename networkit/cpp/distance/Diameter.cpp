// no-networkit-format
/*
 * Diameter.cpp
 *
 *  Created on: 19.02.2014
 *      Author: Daniel Hoske, Christian Staudt
 */

#include <numeric>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/Eccentricity.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

Diameter::Diameter(const Graph& G, DiameterAlgo algo, double error, count nSamples) : Algorithm(), G(&G), error(error), nSamples(nSamples) {
    if (algo == DiameterAlgo::automatic) {
        this->algo = DiameterAlgo::exact;
    } else {
        this->algo = algo;
        if (this->algo == DiameterAlgo::estimatedRange) {
            if (error == -1.f) throw std::invalid_argument("For Diameter::estimatedRange the parameter error(>=0) has to be supplied");
        } else if (this->algo == DiameterAlgo::estimatedSamples) {
            if (nSamples == 0) throw std::invalid_argument("For Diameter::estimatedSamples the parameter nSamples(>0) has to be supplied");
        }
    }
}

void Diameter::run() {
    diameterBounds = {0, 0};
    if (algo == DiameterAlgo::exact) {
        std::get<0>(diameterBounds) = this->exactDiameter(*G);
    } else if (algo == DiameterAlgo::estimatedRange) {
        diameterBounds = this->estimatedDiameterRange(*G, error);
    } else if (algo == DiameterAlgo::estimatedSamples) {
        std::get<0>(diameterBounds) = this->estimatedVertexDiameter(*G, nSamples);
    } else if (algo == DiameterAlgo::estimatedPedantic) {
        std::get<0>(diameterBounds) = this->estimatedVertexDiameterPedantic(*G);
    } else {
        throw std::runtime_error("should never reach this code as the algorithm should be set correctly in the constructor or fail there");
    }
    hasRun = true;
}

std::pair<count, count> Diameter::getDiameter() const {
    return diameterBounds;
}

edgeweight Diameter::exactDiameter(const Graph& G) {
    using namespace std;

    Aux::SignalHandler handler;

    edgeweight diameter = 0.0;

    if (! G.isWeighted()) {
        std::tie(diameter, std::ignore) = estimatedDiameterRange(G, 0);
    } else {
        G.forNodes([&](node v) {
            handler.assureRunning();
            Dijkstra dijkstra(G, v);
            dijkstra.run();
            auto distances = dijkstra.getDistances();
            G.forNodes([&](node u) {
                if (diameter < distances[u]) {
                    diameter = distances[u];
                }
            });
        });
    }

    if (diameter == std::numeric_limits<edgeweight>::max()) {
        throw std::runtime_error("Graph not connected - diameter is infinite");
    }
    return diameter;
}

std::pair<edgeweight, edgeweight> Diameter::estimatedDiameterRange(const Graph &G, double error) {
    if (G.isDirected()) {
        throw std::runtime_error("Error, the diameter of directed graphs cannot be computed yet.");
    }
    if (G.isWeighted()) {
        WARN("The input graph is weighted, but this algorithm ignores weights.");
    }

    Aux::SignalHandler handler;
    /*
     * This is an implementation of a slightly modified version of the exactSumSweep algorithm as presented in
     * Fast diameter and radius BFS-based computation in (weakly connected) real-world graphs: With an application to the six degrees of separation games
     * by Michele Borassi, Pierluigi Crescenzi, Michel Habib, Walter A. Kosters, Andrea Marino, Frank W. Takes
     * http://www.sciencedirect.com/science/article/pii/S0304397515001644
     */

    std::vector<count> eccLowerBound(G.upperNodeIdBound()), eccUpperBound(G.upperNodeIdBound());
    std::vector<bool> finished(G.upperNodeIdBound());

    for (node u = 0; u < G.upperNodeIdBound(); ++u) {
        if (G.hasNode(u)) {
            eccUpperBound[u] = G.numberOfNodes();
        }
    }

    ConnectedComponents comp(G);
    comp.run();
    count numberOfComponents = comp.numberOfComponents();


    std::vector<count> distFirst;
    std::vector<count> ecc(numberOfComponents, 0);
    std::vector<count> distances(G.upperNodeIdBound(), 0);

    count numBFS = 0;

    auto runBFS = [&](const std::vector<node> &startNodes) {
        ++numBFS;
        distFirst.clear();
        distFirst.resize(numberOfComponents, 0);
        std::vector<bool> foundFirstDeg2Node(numberOfComponents, false);

        Traversal::BFSfrom(G, startNodes.begin(), startNodes.end(), [&](node v, count dist) {
            distances[v] = dist;

            index c = comp.componentOfNode(v);
            ecc[c] = std::max(dist, ecc[c]);

            if (!foundFirstDeg2Node[c] && G.degree(v) > 1) {
                foundFirstDeg2Node[c] = true;
                distFirst[c] = dist;
            }
        });

        G.forNodes([&](node u) {
            if (finished[u]) return;

            auto c = comp.componentOfNode(u);

            auto eccValue = std::max(distances[u], ecc[c] - distances[u]);
            eccLowerBound[u] = std::max(eccValue, eccLowerBound[u]);

            if (distances[u] <= distFirst[c]) {
                eccUpperBound[u] = eccValue;
            } else {
                eccUpperBound[u] = std::min(eccUpperBound[u], distances[u] + ecc[c] - 2 * distFirst[c]);
            }

            finished[u] = (eccUpperBound[u] == eccLowerBound[u]);
        });

        ecc.clear();
        ecc.resize(numberOfComponents, 0);
    };

    auto diameterBounds = [&]() {
        count maxExact = *std::max_element(eccLowerBound.begin(), eccLowerBound.end());
        count maxPotential = *std::max_element(eccUpperBound.begin(), eccUpperBound.end());
        return std::make_pair(maxExact, maxPotential);
    };

    std::vector<node> startNodes(numberOfComponents, 0), maxDist(numberOfComponents, 0);
    std::vector<count> maxDeg(numberOfComponents, 0);

    count lb = 0, ub = G.numberOfNodes();

    auto printStartNodes = [&]() {
        DEBUG("Start nodes (lb = ", lb, ", ub = ", ub, "): ");
        for (node u : startNodes) {
            (void)u; // prevent unused variable warning
            DEBUG("Node ", u, " with distance ", distances[u], ", lower bound ", eccLowerBound[u], ", upper bound ", eccUpperBound[u]);
        }
    };

    // for each component, find the node with the maximum degree and add it as start node
    G.forNodes([&](node v) {
        count d = G.degree(v);
        count c = comp.componentOfNode(v);
        if (d >= maxDeg[c]) {
            startNodes[c] = v;
            maxDeg[c] = d;
        }
    });

    handler.assureRunning();

    runBFS(startNodes);

    std::tie(lb, ub) = diameterBounds();

    for (index i = 0; i < 2*G.numberOfNodes() && ub > (lb + error*lb); ++i) {
        handler.assureRunning();
        startNodes.clear();
        startNodes.resize(numberOfComponents, none);

        if ((i % 2) == 0) {
            G.forNodes([&](node u) {
                auto c = comp.componentOfNode(u);
                if (startNodes[c] == none || std::tie(eccUpperBound[u], distances[u]) > std::tie(eccUpperBound[startNodes[c]], distances[startNodes[c]])) {
                    startNodes[c] = u;
                }
            });
        } else {
            G.forNodes([&](node u) {
                auto c = comp.componentOfNode(u);
                // Idea: we select a node that is central (i.e. has a low lower bound) but that is also close to the previous, non-central node.
                // More generally, the best upper bound we can hope for a node v is eccLowerBound[u] + distance(u, v).
                // We select the node the provides the best upper bound for the previous node u in the hope that in its neighborhood there are more nodes for which the bounds can be decreased.
                // Among all these nodes we select the one that has the largest distance to the previous start node.
                if (startNodes[c] == none) {
                    startNodes[c] = u;
                } else {
                    auto compU = eccLowerBound[u] + distances[u], compStart = eccLowerBound[startNodes[c]] + distances[startNodes[c]];
                    if (distances[u] > distFirst[c] && (compU < compStart || (compU == compStart && distances[u] > distances[startNodes[c]]))) {
                        startNodes[c] = u;
                    }
                }
            });
        }

        handler.assureRunning();

        printStartNodes();
        runBFS(startNodes);

        std::tie(lb, ub) = diameterBounds();
    }

    INFO(numBFS, " BFS used");

    return {lb, ub};
}

edgeweight Diameter::estimatedVertexDiameter(const Graph& G, count samples) {

    edgeweight infDist = std::numeric_limits<edgeweight>::max();

    // TODO: consider weights

    auto estimateFrom = [&](node v) -> count {
        BFS bfs(G, v);
        bfs.run();
        auto distances = bfs.getDistances();

        // get two largest path lengths
        edgeweight maxD = 0;
        edgeweight maxD2 = 0; // second largest distance
        for (auto d : distances) {
            if ((d != infDist) && (d >= maxD)) {
                maxD2 = maxD;
                maxD = d;
            }
        }

        edgeweight dMax = maxD + maxD2;
        count vd = (count) dMax + 1; 	// count the nodes, not the edges
        return vd;
    };

    edgeweight vdMax = 0;
    #pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(samples); ++i) {
        node u = GraphTools::randomNode(G);
        edgeweight vd = estimateFrom(u);
        DEBUG("sampled vertex diameter from node ", u, ": ", vd);
        #pragma omp critical
        {
            if (vd > vdMax) {
                vdMax = vd;
            }
        }
    }

    return vdMax;
}

edgeweight Diameter::estimatedVertexDiameterPedantic(const Graph& G) {
    count vd = 0;
    if (!G.isWeighted()) {
        std::vector<bool> visited(G.upperNodeIdBound(), false);
        // perform breadth-first searches
        G.forNodes([&](node u) {
            if (visited[u] == false) {
                count maxDist = 0, maxDist2 = 0;
                Traversal::BFSfrom(G, u, [&](node v, count dist) {
                    visited[v] = true;
                    if (dist > maxDist) {
                        maxDist2 = maxDist;
                        maxDist = dist;
                    }
                    else if (dist > maxDist2) {
                        maxDist2 = dist;
                    }
                });
                if (maxDist + maxDist2 > vd) {
                    vd = maxDist + maxDist2;
                }
                assert (visited[u] == true);
            }
        });
        vd ++; //we need the number of nodes, not the number of edges
    }
    else {
        ConnectedComponents cc(G);
        DEBUG("finding connected components");
        cc.run();
        INFO("Number of components ", cc.numberOfComponents());
        DEBUG("Estimating size of the largest component");
        std::map<count, count> sizes = cc.getComponentSizes();
        count largest_comp_size = 0;
        for(auto it = sizes.cbegin(); it != sizes.cend(); ++it) {
            DEBUG(it->second);
            if (it->second > largest_comp_size) {
                largest_comp_size = it->second;
            }
        }
        INFO("Largest component size: ", largest_comp_size);
        vd = largest_comp_size;
    }
    return vd;
}

std::string Diameter::toString() const {
    return "Diameter()";
}

} /* namespace NetworKit */
