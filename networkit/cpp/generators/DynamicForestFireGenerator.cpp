// no-networkit-format
/*
 * DynamicForestFireGenerator.cpp
 *
 *  Created on: 17.01.2014
 *      Author: cls
 */

#include <queue>
#include <set>
#include <unordered_map>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/generators/DynamicForestFireGenerator.hpp>

namespace NetworKit {

using neighborFunction = std::function<std::vector<node>(node)>;

DynamicForestFireGenerator::DynamicForestFireGenerator(double p, bool directed, double r) :p(p), directed(directed), r(r), firstCall(true) {
    G = Graph(0, false, directed);
}

std::vector<GraphEvent> DynamicForestFireGenerator::generate(count nSteps) {

    std::vector<GraphEvent> stream;
    std::set<node> empty;

    /* this function creates a new node and connects it to
     * other nodes according to the forest fire model
     */
    auto connectNewNode = [&]() {
        //list of nodes that were visited
        std::unordered_map<node, bool> visited;
        // nodes which were found but not processed
        std::queue<node> activeNodes;
        /* vector of nodes that visited
         * and the new node will connect to
         */
        std::vector<node> burnedNodes;

        auto forwardNeighbors = [&](node u) {
            std::vector<node> validEdges;
            G.forNeighborsOf(u, [&](node x){
                if (! visited[x]) {
                    validEdges.push_back(x);
                }
            });
            return validEdges;
        };

        auto backwardNeighbors = [&](node u) {
            std::vector<node> validEdges;
            G.forInNeighborsOf(u, [&](node, node x){
                if (! visited[x]) {
                    validEdges.push_back(x);
                }
            });
            return validEdges;
        };

        // select "edges" of node u
        auto selectEdges = [&](node u, double prob, const neighborFunction &getNeighbors) {
            /* combine all valid edges (edges to non-visited nodes)
             * into a vector that we can randomly select one
             */
            std::vector<node> validEdges = getNeighbors(u);
            std::set<node> edges;
            while (true) {
                /* get geometric distribution by burning edges
                 * until first failure
                 */
                double q = Aux::Random::real(1.0);
                if (q > prob || validEdges.empty()) {
                    break;
                }
                count index = Aux::Random::integer(validEdges.size() - 1);
                edges.insert(validEdges[index]);
                validEdges[index] = validEdges.back();
                validEdges.pop_back();
            }
            return edges;
        };

        // select ambassador node
        node a = none;
        do {
            a = Aux::Random::integer(G.upperNodeIdBound());
        } while (! G.hasNode(a));
        assert (a != none);
        DEBUG("selected ambassador: ", a);

        node v = G.addNode();
        stream.emplace_back(GraphEvent::NODE_ADDITION, v);
        DEBUG("created node ", v);

        visited[a] = true;
        activeNodes.push(a);
        burnedNodes.push_back(a);

        // burn through the graph in a BFS-like fashion
        while (! activeNodes.empty()) {
            node w = activeNodes.front();
            activeNodes.pop();
            std::set<node> edges = selectEdges(w, p, forwardNeighbors);;
            if (directed) {
                std::set<node> backwardEdges = selectEdges(w, p*r, backwardNeighbors);
                edges.insert(backwardEdges.begin(), backwardEdges.end());
            }
            for (node x : edges) {
                activeNodes.push(x);
                burnedNodes.push_back(x);
                visited[x] = true;
            }
        }

        for (node w : burnedNodes) {
            G.addEdge(v, w);
            stream.emplace_back(GraphEvent::EDGE_ADDITION, v, w);
        }
    };

    // initial graph
    if (firstCall && nSteps > 0) {
        node s = G.addNode();
        stream.emplace_back(GraphEvent::NODE_ADDITION, s);
        stream.emplace_back(GraphEvent::TIME_STEP);
        firstCall = false;
        --nSteps;
    }

    for (index step = 0; step < nSteps; ++step) {
        connectNewNode();
        stream.emplace_back(GraphEvent::TIME_STEP);
    }

    return stream;
}

} /* namespace NetworKit */
