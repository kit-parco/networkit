// no-networkit-format
/*
 * DynBetweenness.cpp
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
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/SSSP.hpp>
#include <networkit/centrality/DynBetweenness.hpp>

namespace NetworKit {

DynBetweenness::DynBetweenness(const Graph &G) : Centrality(G),
distances(G.upperNodeIdBound(), std::vector<edgeweight>(G.upperNodeIdBound())),
distancesOld(G.upperNodeIdBound(), std::vector<edgeweight>(G.upperNodeIdBound())),
sigma(G.upperNodeIdBound(), std::vector<edgeweight>(G.upperNodeIdBound())),
sigmaOld(G.upperNodeIdBound(), std::vector<edgeweight>(G.upperNodeIdBound())){}

/**
 * Run method that stores a single shortest path for each node pair and stores shortest distances
 */
void DynBetweenness::run() {
    count z = G.upperNodeIdBound();
    scoreData.clear();
    scoreData.resize(z);

    auto computeDependencies = [&](node s) {
            // run SSSP algorithm and keep track of everything
            std::unique_ptr<SSSP> sssp;
            if (G.isWeighted()) {
                    sssp = std::make_unique<Dijkstra>(G, s, true, true);
            } else {
                    sssp = std::make_unique<BFS>(G, s, true, true);
            }

            sssp->run();

            G.forNodes([&](node t){
                    distances[s][t] = sssp->distance(t);
                    sssp->numberOfPaths(t).ToDouble(sigma[s][t]);
            });

            // compute dependencies for nodes in order of decreasing distance from s
            std::vector<node> stack = sssp->getNodesSortedByDistance();
            // set maxDistance to the distance of the furthest vertex
            if(distances[s][stack.back()] > diameter){
                diameter = distances[s][stack.back()];
            }
            std::vector<double> dependencies(z, 0);
            while (!stack.empty()) {
                    node t = stack.back();
                    stack.pop_back();
                    G.forInNeighborsOf(t, [&](node p, edgeweight edgept) {
                            if (distances[s][t] == distances[s][p] + edgept) {
                                double weight = sigma[s][p] / sigma[s][t];
                                dependencies[p] += weight * (1 + dependencies[t]);
                            }
                    });
                    TRACE("Dependency of node ",s, " on node ", t, ": ", dependencies[t]);
                    if (t != s) {
                            scoreData[t] += dependencies[t];
                    }
            }
    };

    G.forNodes(computeDependencies);
    distancesOld = distances;
    sigmaOld = sigma;
    hasRun = true;
}


void DynBetweenness::increaseScore(std::vector<bool> & affected, node y, std::priority_queue<std::pair<double, node>, std::vector<std::pair<double,node>>, CompareDist> & Q) {
    std::vector<double> dep(G.upperNodeIdBound(), 0);
    std::vector<bool> visited(G.upperNodeIdBound(), false);
    while (!Q.empty()) {
        affectedDep ++;
    //	TRACE("(Before )Size: ", Q.size());
        node x = Q.top().second; // notice that the keys are diam - distance, so we actually extract in order of decreasing distance
        Q.pop();
        TRACE("Extracted node ", x);
    //	TRACE("(After )Size: ", Q.size());
        scoreData[x] += dep[x];
        TRACE("Dependency of ",y, " on ",x,": ", dep[x]);
        if (!G.isDirected()) {
            scoreData[x] += dep[x];
        }
        G.forNeighborsOf(x, [&](node w, edgeweight weightxw){
            if (w!= x && x != y && distances[x][y] == distances[w][y] + weightxw) {
                if (affected[x]) {
                    TRACE("Affected node ", x, ". Adding to dep of ", w, ": ", sigma[w][y]/sigma[x][y]*(1+dep[x]));
                    dep[w] += sigma[w][y]/sigma[x][y]*(1+dep[x]);
                } else {
                    TRACE("Non affected node ", x, ". Adding to dep of ", w, ": ", sigma[w][y]/sigma[x][y]*(dep[x]));
                    dep[w] += sigma[w][y]/sigma[x][y]*(dep[x]);
                }
                if (! visited[w] && !affected[w] && w != y) {
                    TRACE("Inserting node ", w, " with new priority ", distances[w][y]);
                    Q.push(std::make_pair(diameter + 1.0 - distances[w][y], w));
                    visited[w] = true;
                }
            }
        });
    }
}


void DynBetweenness::decreaseScore(std::vector<bool> & affected, node y, std::priority_queue<std::pair<double, node>, std::vector<std::pair<double,node>>, CompareDist> & Q) {
    std::vector<double> dep(G.upperNodeIdBound(), 0);
    std::vector<bool> visited(G.upperNodeIdBound(), false);
    while (!Q.empty()) {
        affectedDep ++;
    //	TRACE("(Before )Size: ", Q.size());
        node x = Q.top().second; // notice that the keys are diam - distance, so we actually extract in order of decreasing distance
        Q.pop();
        TRACE("Extracted node ", x);
    //	TRACE("(After )Size: ", Q.size());
        scoreData[x] -= dep[x];
        // TRACE("Dependency of ",y, " on ",x,": ", dep[x]);
        if (!G.isDirected()) {
            scoreData[x] -= dep[x];
        }
        G.forNeighborsOf(x, [&](node w, edgeweight weightxw){
            if (w!= x && (w != v || x != u) && x != y && distancesOld[x][y] == distancesOld[w][y] + weightxw) {
                if(sigmaOld[x][y] > 0){
                    if (affected[x]) {
                        // TRACE("Affected node ", x, ". Subtracting from dep of ", w, ": ", sigmaOld[w][y]/sigmaOld[x][y]*(1+dep[x]));
                        dep[w] += sigmaOld[w][y]/sigmaOld[x][y]*(1+dep[x]);
                    } else {
                        // TRACE("Non affected node ", x, ". Subtracting from dep of ", w, ": ", sigmaOld[w][y]/sigmaOld[x][y]*(dep[x]));
                        dep[w] += sigmaOld[w][y]/sigmaOld[x][y]*(dep[x]);
                    }
                }
                if (! visited[w] && !affected[w] && w != y) {
                    TRACE("Inserting node ", w, " with old priority ", distancesOld[w][y]);
                    Q.push(std::make_pair(diameter + 1.0 - distancesOld[w][y], w));
                    visited[w] = true;
                }
            }
        });
    }
}


void DynBetweenness::update(GraphEvent event) {
    timeDep = 0;
    INFO("Diameter: ", diameter);
    visitedPairs = 0;
    INFO("Entering update");
    u = event.u;
    v = event.v;
    edgeweight weightuv = G.weight(u,v);
    if (!(event.type==GraphEvent::EDGE_ADDITION || (event.type==GraphEvent::EDGE_WEIGHT_INCREMENT && event.w < 0))) {
        throw std::runtime_error("event type not allowed. Edge insertions and edge weight decreases only.");
    }
    if (weightuv < distances[u][v]) {
        // initializations
        affectedAPSP = 0;
        affectedDep = 0;
        INFO("Old distance: ", distances[u][v]);
        std::queue<std::pair<node, node>> modified;
        count z = G.upperNodeIdBound();
        std::vector<node> source_nodes(z);
        std::vector<node> n_sources(z, 0);
        std::queue<node> Q;
        // phase 1: find affected source nodes using bfs
        count i = 0;
        std::queue<node> bfsQ;
        std::vector<bool> visited(z, false);
        INFO("Phase 1. distances[", u,"][", v,"] = ", distances[u][v], ", and G.weight", u,", ", v," = ",G.weight(u,v));
        distances[u][v] = weightuv;
        modified.push(std::make_pair(u, v));
        sigma[u][v] = 1;
        visited[u] = true;
        if(!G.isDirected()) {
            distances[v][u] = distances[u][v];
            sigma[v][u] = 1;
        }
        bfsQ.push(u);
        INFO("Entering bfs");
        while (! bfsQ.empty()) {
            node x = bfsQ.front();
            bfsQ.pop();
            DEBUG("Dequeueing node ", x);
            G.forInNeighborsOf(x, [&](node w, edgeweight) { // identify and process neighbors w of x
                if (visited[w] == false && distances[w][v] >= distances[w][u] + weightuv) {
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
        std::vector<bool> enqueued(G.upperNodeIdBound(), false);
        enqueued[v] = true;

        while (! stack.empty()) {
            node y = stack.top();
            if (!visited[y]) {
                std::priority_queue<std::pair<double, node>, std::vector<std::pair<double,node>>, CompareDist> Qnew;
                std::priority_queue<std::pair<double, node>, std::vector<std::pair<double,node>>, CompareDist> Qold;
                std::vector<bool> affected(G.upperNodeIdBound(), false);
                affected[u] = true;
                // we leave y in the stack (so that we know when we're done visiting the subtree rooted in y)
                n_sources[y] = n_sources[Pred[y]];
                visited[y] = true;
                // since u is not in source, we insert it now
                // Qnew.insert(diameter + 1 - distances[u][y], u);
                // Qold.insert(diameter + 1 - distancesOld[u][y], u);
                Qnew.push(std::make_pair(diameter + 1.0 - distances[u][y], u));
                Qold.push(std::make_pair(diameter + 1.0 - distancesOld[u][y], u));
                for (count c = 0; c < n_sources[y]; c++) {
                    node s = source_nodes[c];
                    if (s != u) {
                        affectedAPSP ++;
                    }
                    if (distances[s][y] > distances[s][u] + weightuv + distances[v][y]) {
                        distances[s][y] = distances[s][u] + weightuv + distances[v][y];
                        sigma[s][y] = sigma[s][u] * sigma[v][y];
                        if(!G.isDirected()) {
                            distances[y][s] = distances[s][y];
                            sigma[y][s] = sigma[s][y];
                        }
                        modified.push(std::make_pair(s, y));
                        affected[s] = true;
                        TRACE("Node ", y,", Inserting node ", s, " with new priority ", diameter + 1 - distances[s][y]);
                        TRACE("Node ", y,", Inserting node ", s, " with old priority ", diameter + 1 - distancesOld[s][y]);
                        // Qnew.insert(diameter + 1 - distances[s][y], s);
                        // Qold.insert(diameter + 1 - distancesOld[s][y], s);
                        Qnew.push(std::make_pair(diameter + 1 - distances[s][y], s));
                        Qold.push(std::make_pair(diameter + 1 - distancesOld[s][y], s));
                    } else if (distances[s][y] == distances[s][u] + weightuv + distances[v][y]) {
                            sigma[s][y] += sigma[s][u] * sigma[v][y];
                            if(!G.isDirected()) {
                                sigma[y][s] = sigma[s][y];
                            }
                            modified.push(std::make_pair(s, y));
                            affected[s] = true;
                            TRACE("Node ", y,", Inserting node ", s, " with new priority ", diameter + 1 - distances[s][y]);
                            TRACE("Node ", y,", Inserting node ", s, " with old priority ", diameter + 1 - distancesOld[s][y]);
                            // Qnew.insert(diameter + 1 - distances[s][y], s);
                            // Qold.insert(diameter + 1 - distancesOld[s][y], s);
                            Qnew.push(std::make_pair(diameter + 1.0 - distances[s][y], s));
                            Qold.push(std::make_pair(diameter + 1.0 - distancesOld[s][y], s));
                    }
                    else if (distances[s][y] < distances[s][u] + weightuv + distances[v][y]) {
                        std::swap(source_nodes[c], source_nodes[n_sources[y] - 1]);
                        c --;
                        n_sources[y] --;
                    }
                }
                // now we update the bc scores
                TRACE("Size of Q old: ", Qold.size());
                TRACE("Size of Q new: ", Qnew.size());
                clock_t tStart = clock();
                increaseScore(affected, y, Qnew);
                decreaseScore(affected, y, Qold);
                timeDep += (double)(clock() - tStart)/CLOCKS_PER_SEC;
                // adding successors of y to the stack
                G.forNeighborsOf(y, [&](node w, edgeweight weightyw){
                    // we go down the BFS tree rooted in v in a DFS order (the last check is necessary to make sure that (y, w) is an edge of the BFS tree rooted in v)
                    if (w!=y && visited[w] == false && enqueued[w] == false && distances[u][w] >= distances[v][w] + weightuv && distances[v][w] == distances[v][y] + weightyw) {
                        if (distances[u][w] > distances[v][w] + weightuv) {
                            distances[u][w] = distances[v][w] + weightuv;
                            TRACE(" > Setting sigma ", u, ",", w, " from  ", sigma[u][w], " to ", sigma[v][w]);
                            sigma[u][w] = sigma[v][w];
                        } else if (distances[u][w] == distances[v][w] + weightuv) {
                            TRACE(" = Increasing sigma ", u, ",", w, " from  ", sigma[u][w], " of ", sigma[v][w]);
                            sigma[u][w] += sigma[v][w];
                        }
                        if(!G.isDirected()) {
                            distances[w][u] = distances[u][w];
                            sigma[w][u] = sigma[u][w];
                        }
                        stack.push(w);
                        enqueued[w] = true;
                        modified.push(std::make_pair(u, w));
                        Pred[w] = y;
                    }
                });
            } else {
                // we remove y from the stack
                stack.pop();
            }
        }

        // reset sigma old to sigma new and distance old to distance new
        while(!modified.empty()) {
            std::pair<node, node> p = modified.front();
            modified.pop();
            node u = p.first;
            node v = p.second;
            // set also the "old" data structures to the new values
            distancesOld[u][v] = distances[u][v];
            sigmaOld[u][v] = sigma[u][v];
            if(!G.isDirected()) {
                distancesOld[v][u] = distances[u][v];
                sigmaOld[v][u] = sigma[u][v];
            }
        }
        // distancesOld = distances;
        // sigmaOld = sigma;

    }

}

void DynBetweenness::updateBatch(const std::vector<GraphEvent>& batch) {
  for(auto e : batch){
    update(e);
  }
}

edgeweight DynBetweenness::getDistance(node u, node v) {
    return distances[u][v];
}

edgeweight DynBetweenness::getSigma(node u, node v) {
    return sigma[u][v];
}

count DynBetweenness::visPairs() {
    return visitedPairs;
}

count DynBetweenness::numAffectedAPSP(){
    return affectedAPSP;
}

count DynBetweenness::numAffectedDep(){
    return affectedDep;
}

double DynBetweenness::getTimeDep(){
    return timeDep;
}


} /* namespace NetworKit */
