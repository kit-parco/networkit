// no-networkit-format
/*
 * DynAPSP.cpp
 *
 *  Created on: 12.08.2015
 *      Author: Arie Slobbe, Elisabetta Bergamini
 */

#include <algorithm>
#include <ctime>
#include <memory>
#include <queue>
#include <unordered_set>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/distance/APSP.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/DynAPSP.hpp>
#include <networkit/distance/SSSP.hpp>

namespace NetworKit {

DynAPSP::DynAPSP(Graph& G) : APSP(G) {}

/**
 * Run method that stores a single shortest path for each node pair and stores shortest distances
 */
void DynAPSP::run() {
    distances.resize(G.upperNodeIdBound());
    G.parallelForNodes([&](node u){
        std::unique_ptr<SSSP> sssp;
        if (G.isWeighted()) {
            sssp = std::make_unique<Dijkstra>(G, u);
        } else {
            sssp = std::make_unique<BFS>(G, u);
        }
        sssp->run();
        distances[u] = sssp->getDistances();
    });
  hasRun = true;
}

std::vector<node> DynAPSP::getPath(node u, node v) {
    std::vector<node> path = {};
    if (distances[u][v] < std::numeric_limits<edgeweight>::max()) {
        node current = v;
        while (current != u) {
            path.push_back(current);
            G.forInEdgesOf(current, [&](node z, edgeweight w) {
                if (distances[u][current] == distances[u][z] + w) {
                    current = z;
                }
            });
        }
        path.push_back(u);
        std::reverse(path.begin(), path.end());
    }
    return path;
}



void DynAPSP::update(GraphEvent event) {
    visitedPairs = 0;
    INFO("Entering update");
    node u = event.u;
    node v = event.v;
    edgeweight weightuv = G.weight(u,v);
    if (!(event.type==GraphEvent::EDGE_ADDITION || (event.type==GraphEvent::EDGE_WEIGHT_INCREMENT && event.w < 0))) {
        throw std::runtime_error("event type not allowed. Edge insertions and edge weight decreases only.");
    }
    if (weightuv < distances[u][v]) {
        // initializations
        count z = G.upperNodeIdBound();
        std::vector<node> source_nodes(z);
        std::vector<node> n_sources(z, 0);
        std::queue<node> Q;
        std::vector<bool> enqueued(G.upperNodeIdBound(), false);
        // phase 1: find affected source nodes using bfs
        count i = 0;
        std::queue<node> bfsQ;
        std::vector<bool> visited(z, false);
        INFO("Phase 1. distances[", u,"][", v,"] = ", distances[u][v], ", and G.weight", u,", ", v," = ",G.weight(u,v));
        distances[u][v] = weightuv;
        if(!G.isDirected()) {
            distances[v][u] = distances[u][v];
        }
        bfsQ.push(u);
        INFO("Entering bfs");
        while (! bfsQ.empty()) {
            node x = bfsQ.front();
            bfsQ.pop();
            DEBUG("Dequeueing node ", x);
            G.forInNeighborsOf(x, [&](node w, edgeweight) { // identify and process neighbors w of x
                if (visited[w] == false && distances[w][v] > distances[w][u] + weightuv) {
                    bfsQ.push(w);
                    DEBUG("Pushing neighbor ", w);
                    visited[w] = true;
                    source_nodes[i] = w;
                    i ++;
                }
            });
        }
        // notice that source nodes does not contain u
        n_sources[u] = i;
        // phase 2: for all source nodes, update distances to affected sinks
        std::vector<node> Pred(G.upperNodeIdBound());
        Pred[v] = u;
        std::stack<node> stack;
        stack.push(v);
        visited.clear();
        visited.resize(z, false);
        while (! stack.empty()) {
            node y = stack.top();
            if (!visited[y]) {
                // we leave y in the stack (so that we know when we're done visiting the subtree rooted in y)
                n_sources[y] = n_sources[Pred[y]];
                visited[y] = true;
                for (count c = 0; c < n_sources[y]; c++) {
                    node s = source_nodes[c];
                    if (distances[s][y] > distances[s][u] + weightuv + distances[v][y]) {
                        distances[s][y] = distances[s][u] + weightuv + distances[v][y];
                        if(!G.isDirected()) {
                            distances[y][s] = distances[s][y];
                        }
                    } else {
                        std::swap(source_nodes[c], source_nodes[n_sources[y] - 1]);
                        c --;
                        n_sources[y] --;
                    }
                }
                // adding successors of y to the stack
                G.forNeighborsOf(y, [&](node w, edgeweight weightyw){
                    // we go down the BFS tree rooted in v in a DFS order (the last check is necessary to make sure that (y, w) is an edge of the BFS tree rooted in v)
                    if (visited[w] == false && distances[u][w] > distances[v][w] + weightuv && distances[v][w] == distances[v][y] + weightyw) {
                        distances[u][w] = distances[v][w] + weightuv;
                        if(!G.isDirected()) {
                            distances[w][u] = distances[u][w];
                        }
                        stack.push(w);
                        Pred[w] = y;
                    }
                });
            } else {
                // we remove y from the stack
                stack.pop();
            }
        }
        // while(!Q.empty()) {
        // 	node y = getMin();
        // 	enqueued[y] = false;
        // 	// update for all source nodes
        // 	for (node x: source_nodes[Pred[y]]) {
        // 		visitedPairs ++;
        // 		if (distances[x][y] > distances[x][u] + distances[u][y]) {
        // 			distances[x][y] = distances[x][u] + distances[u][y];
        // 			if(!G.isDirected()) {
        // 				distances[y][x] = distances[x][y];
        // 			}
        // 			source_nodes[y].push_back(x);
        // 		}
        // 	}
        //
        // 	G.forNeighborsOf(y, [&](node w, edgeweight weightyw){
        // 		if (distances[u][w] > distances[u][y] + weightyw && distances[v][w] == distances[v][y] + weightyw) { // I also check that y was a predecessor for w in the s.p. from v
        // 			distances[u][w] = distances[u][y] + weightyw;
        // 			Pred[w] = y;
        // 			if(!G.isDirected()) {
        // 				distances[w][u] = distances[u][w];
        // 			}
        // 			updateQueue(w, distances[u][w]);
        // 		}
        // 	});
        // }
    }
}

void DynAPSP::updateBatch(const std::vector<GraphEvent>& batch) {
  for(auto e : batch){
    update(e);
  }
}

count DynAPSP::visPairs() {
    return visitedPairs;
}


} /* namespace NetworKit */
