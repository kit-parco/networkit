// no-networkit-format
/*
 * DynBetweennessOneNode.cpp
 *
 *  Created on: 10.03.2016
 *      Author: Elisabetta Bergamini
 */

#include <algorithm>
#include <queue>
#include <unordered_set>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/PrioQueue.hpp>
#include <networkit/centrality/DynBetweennessOneNode.hpp>

namespace NetworKit {

DynBetweennessOneNode::DynBetweennessOneNode(Graph& G, node x) : G(G), x(x) {}

/**
 * Run method that stores a One shortest path for each node pair and stores shortest distances
 */
void DynBetweennessOneNode::run() {
    bcx = 0;
    distances.resize(G.upperNodeIdBound());
    sigma.resize(G.upperNodeIdBound());
    sigmax.resize(G.upperNodeIdBound());
    Pred.resize(G.upperNodeIdBound());
  // for each node, we execute a BFS and compute number of SPs between node pairs and number of SPs between node pairs that go through x
    G.forNodes([&](node source){
        distances[source].resize(G.upperNodeIdBound(), infDist);
        sigma[source].resize(G.upperNodeIdBound(), 0);
        sigmax[source].resize(G.upperNodeIdBound(), 0);
        std::queue<node> q;
        std::vector<bool> visited(G.upperNodeIdBound(), false);
        q.push(source);
        visited[source] = true;
        sigma[source][source] = 1;
        distances[source][source] = 0;
        while (! q.empty()) {
            node u = q.front();
            q.pop();
            if (u == x) {
                sigmax[source][u] = sigma[source][u];
            }
            // insert untouched neighbors into queue
            G.forNeighborsOf(u, [&](node v) {
                if (!visited[v]) {
                    q.push(v);
                    visited[v] = true;
                    distances[source][v] = distances[source][u] + 1;
                }
                if (distances[source][v] == distances[source][u] + 1) {
                    sigma[source][v] += sigma[source][u]; 	// all the shortest paths to u are also shortest paths to v now
                    if (u == x) {
                        sigmax[source][v] = sigma[source][u];
                    } else {
                        sigmax[source][v] += sigmax[source][u]; // it is 0 until we visit x
                    }
                }
            });
        }
    });

    G.forNodes([&](node s){
        G.forNodes([&](node t){
            if (t != x && s != x && sigma[s][t] != 0) {
                bcx += double(sigmax[s][t])/sigma[s][t];
            }
        });
    });
}



void DynBetweennessOneNode::update(GraphEvent event) {
    // DEBUG("Entering update");
    node u = event.u;
    node v = event.v;
    edgeweight weightuv = G.weight(u,v);
    if (!(event.type==GraphEvent::EDGE_ADDITION || (event.type==GraphEvent::EDGE_WEIGHT_INCREMENT && event.w < 0))) {
        throw std::runtime_error("event type not allowed. Edge insertions and edge weight decreases only.");
    }
    if (weightuv <= distances[u][v]) {
        // initializations
        count z = G.upperNodeIdBound();
        std::vector<std::vector<node>> source_nodes(z);
        std::queue<node> Q;
        std::vector<bool> enqueued(G.upperNodeIdBound(), false);
        // queue with all visited nodes
        // if u has a new shortest path going through v, it updates the distance of u
        // and inserts u in the priority queue (or updates its priority, if already in Q)
        auto updateQueue = [&](node u) {
            if (!enqueued[u]) {
                Q.push(u);
                enqueued[u] = true;
            }
        };
        // returns smallest element in Q
        auto getMin = [&]() {
            node s = Q.front();
            Q.pop();
            return s;
        };
        // phase 1: find affected source nodes using bfs
        std::queue<node> bfsQ;
        std::vector<bool> visited(z, false);
        INFO("Phase 1. distances[", u,"][", v,"] = ", distances[u][v], ", and G.weight", u,", ", v," = ",G.weight(u,v));
    //	distances[u][v] = weightuv;
        // if(!G.isDirected()) {
        // 	distances[v][u] = distances[u][v];
        // }
        bfsQ.push(u);
        source_nodes[u].push_back(u);
        visited[u] = true;
        INFO("Entering bfs");
        while (! bfsQ.empty()) {
            node s = bfsQ.front();
            bfsQ.pop();
            DEBUG("Dequeueing node ", s);
            G.forInNeighborsOf(s, [&](node w, edgeweight) { // identify and process neighbors w of s
                if (visited[w] == false && distances[w][v] >= distances[w][u] + weightuv) {
                    bfsQ.push(w);
                    DEBUG("Pushing neighbor ", w);
                    visited[w] = true;
                    source_nodes[u].push_back(w);
                }
            });
        }
        // phase 2: for all source nodes, update distances to affected sinks
        std::vector<node> Pred(G.upperNodeIdBound());
        Pred[v] = u;
        visited[u] = true; // what's this?
        updateQueue(v);
        INFO("Entering phase 2. source_nodes[", u,"] = ", source_nodes[u]);
        while(!Q.empty()) {
            node y = getMin();
            enqueued[y] = false;
            // update for all source nodes
            for (node s: source_nodes[Pred[y]]) {
                if (distances[s][y] >= distances[s][u] + weightuv + distances[v][y]) {
                    if (s != x && y!= x && sigma[s][y]) {
                        bcx -= double(sigmax[s][y])/sigma[s][y];
                        if(!G.isDirected()) {
                            bcx -= double(sigmax[s][y])/sigma[s][y];
                        }
                    }
                    if (s == u && y == v) {
                        sigma[u][v] = 1;
                        if (s == x || y == x) {
                            sigmax[u][v] = 1; // TODO fix for weighted graphs, the new edge might be of equal length
                        } else {
                            sigmax[u][v] = 0;
                        }
                    }
                    else {
                        if (distances[s][y] > distances[s][u] + weightuv + distances[v][y]) {
                            if (s ==1 && y ==2 ) INFO(" GREATER s = ", s, ", y = ", y);
                            sigma[s][y] = sigma[s][u]*sigma[v][y];
                            sigmax[s][y] = sigmax[s][u]*sigma[v][y] + sigma[s][u]*sigmax[v][y];
                        } else {
                            if (s ==1 && y ==2 ) INFO(" EQUAL s = ", s, ", y = ", y);
                            sigma[s][y] += sigma[s][u]*sigma[v][y];
                            sigmax[s][y] += sigmax[s][u]*sigma[v][y] + sigma[s][u]*sigmax[v][y];
                        }
                    }
                    if (s != x && y!= x && sigma[s][y]) {
                        bcx += double(sigmax[s][y])/sigma[s][y];
                        if(!G.isDirected()) {
                            bcx += double(sigmax[s][y])/sigma[s][y];
                        }
                    }
                    distances[s][y] = distances[s][u] + weightuv + distances[v][y];
                    if(!G.isDirected()) {
                        distances[y][s] = distances[s][y];
                        sigma[y][s] = sigma[s][y];
                        sigmax[y][s] = sigmax[s][y];
                    }
                    source_nodes[y].push_back(s);
                }
            }
            // loop over all neighbors
            G.forNeighborsOf(y, [&](node w, edgeweight weightyw){
                if (distances[u][w] >= weightuv + distances[v][y] + weightyw && distances[v][w] == distances[v][y] + weightyw) { // I also check that y was a predecessor for w in the s.p. from v
                    Pred[w] = y;
                    updateQueue(w);
                }
            });
        }
    }
}

void DynBetweennessOneNode::updateBatch(const std::vector<GraphEvent>& batch) {
  for(auto e : batch){
    update(e);
  }
}

edgeweight DynBetweennessOneNode::getDistance(node u, node v) {
    return distances[u][v];
}

edgeweight DynBetweennessOneNode::getSigma(node u, node v) {
    return sigma[u][v];
}

edgeweight DynBetweennessOneNode::getSigmax(node u, node v) {
    return sigmax[u][v];
}

edgeweight DynBetweennessOneNode::getbcx() {
    return bcx;
}

} /* namespace NetworKit */
