/*
 * DynApproxBetweenness.cpp
 *
 *  Created on: 31.07.2014
 *      Author: ebergamini
 */

#include "DynApproxBetweenness2.h"
#include "../auxiliary/Random.h"
#include "../properties/Diameter.h"
#include "../graph/Sampling.h"
#include "../graph/DynDijkstra.h"
#include "../graph/DynBFS.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/NumericTools.h"
#include <queue>
#include <math.h>

namespace NetworKit {

DynApproxBetweenness2::DynApproxBetweenness2(const Graph& G, double epsilon, double delta, bool storePredecessors) : Centrality(G, true),
storePreds(storePredecessors), epsilon(epsilon), delta(delta) {

}


count DynApproxBetweenness2::getNumberOfSamples() {
    return r;
}


void DynApproxBetweenness2::run() {
    if (G.isDirected()) {
        throw std::runtime_error("Invalid argument: G must be undirected.");
    }
    scoreData.clear();
    scoreData.resize(G.upperNodeIdBound());
    u.clear();
    v.clear();
    sampledPaths.clear();
    vis.clear();
    vis.resize(G.upperNodeIdBound());
    npaths.clear();
    maxDist.clear();
    maxDist2.clear();


    vd = Diameter::estimatedVertexDiameterPedantic(G);

    INFO("estimated diameter: ", vd);
    r = ceil((c / (epsilon * epsilon)) * (floor(log2(vd - 2)) + 1 + log(1 / delta)));
    INFO("taking ", r, " path samples");
    sssp.clear();

    sampleNewPaths(0, r);
    //we keep track of additional sssp for the components that have never been visited
    std::queue<node> U;
    G.forNodes([&](node n) {
        if (vis[n] == 0) {
            U.push(n);
        }
    });
    sssp2.clear();
    r2 = 0;
    while(!U.empty()) {
        node s = U.front();
        U.pop();
        if (vis[s] > 0) {
            continue;
        }
        BFSvisit bfs(G,s);
        count max1 = 0, max2 = 0;
        bfs.run([&](node w, count dist){
            vis[w]++;
            if (dist > max1) {
                max2 = max1;
                max1 = dist;
            } else if (dist > max2) {
                max2 = dist;
            }
        },
        [&](node w, node z, count dist_w, count dist_z, bool visited){
        });
        sssp2.push_back(bfs);
        maxDist.push_back(max1);
        maxDist2.push_back(max2);
        r2 ++;
    }
    assert(sampledPaths.size() == r);
    assert(npaths.size() == r);
    assert(maxDist.size() == r+r2);
    assert(maxDist2.size() == r+r2);
    hasRun = true;
}


void DynApproxBetweenness2::sampleNewPaths(count start, count end) {
    for (count i = start; i < end; i++) {
        DEBUG("sample ", i);
        // sample random node pair
        node u1, u2;
        u1 = Sampling::randomNode(G);
        do {
            u2 = Sampling::randomNode(G);
        } while (u1 == u2);
        u.push_back(u1);
        v.push_back(u2);
        if (G.isWeighted()) {
            return; //TODO weighted graphs!
        //    sssp[i].reset(new DynDijkstra(G, u[i], storePreds));
        } else {
            BFSvisit bfs(G, u[i]);
            sssp.push_back(bfs);
        }
        DEBUG("running shortest path algorithm for node ", u[i]);
        std::vector<count> sigma(G.upperNodeIdBound(), 0);
        sigma[u[i]] = 1;
        count max1 = 0, max2 = 0;
        sssp[i].run([&](node w, count dist) {
            vis[w]++;
            if (dist > max1) {
                max2 = max1;
                max1 = dist;
            } else if (dist > max2) {
                max2 = dist;
            }
        },
        //we build sigma
        [&](node w, node z, count dist_w, count dist_z, bool visited) {
            if (!visited) {
                sigma[z] = sigma[w];
            } else if (dist_z == dist_w + 1) {
                sigma[z] += sigma[w];
            }
        });
        npaths.push_back(sigma);
        maxDist.push_back(max1);
        maxDist2.push_back(max2);
        std::vector<node> path;
        if (sssp[i].distance(v[i]) < infDist) { // at least one path between {u, v} exists
            // random path sampling and estimation update
            node t = v[i];
            while (t != u[i])  {
                // sample z in P_u(t) with probability sigma_uz / sigma_us
                std::vector<std::pair<node, double>>choices;
                G.forEdgesOf(t, [&](node t, node z, edgeweight w){
                    if (Aux::NumericTools::logically_equal(sssp[i].distance(t), sssp[i].distance(z) + w)) {
                        // workaround for integer overflow in large graphs
                        double weight = npaths[i][z] / double(npaths[i][t]);
                        assert(weight != 0);
                        choices.emplace_back(z, weight);
                    }

                });
                assert (choices.size() > 0);
                node z = Aux::Random::weightedChoice(choices);
                assert (z <= G.upperNodeIdBound());
                if (z != u[i]) {
                    scoreData[z] += 1 / (double) r;
                    path.push_back(z);
                }
                t = z;
            }
        }
        sampledPaths.push_back(path);
    }
}

void DynApproxBetweenness2::update(const std::vector<GraphEvent>& batch) {
    std::queue<node> U;
    count new_vd = vd;
    INFO("Updating the r paths");
    assert(sssp.size() == r);
    for (node i = 0; i < r; i++) {
        DEBUG("Updating distance of source node number ",i," out of ",r);
        assert(npaths.size() == r);
        assert(vis.size() == G.upperNodeIdBound());
        assert(maxDist.size() >= r);
        assert(maxDist2.size() >= r);
        sssp[i].update(batch,
        [&](node w, count dist, count m){
            DEBUG("Updating node ",w,", source node: ",u[i], ", distance: ", dist, ", m: ", m);
            assert(npaths[i].size() == G.upperNodeIdBound());
            npaths[i][w] = 0;
            if(dist == infDist) { //node w has become reachable
                vis[w]++;
            }
            //calculate the maxDist
            if(m > maxDist[i]) {
                maxDist2[i] = maxDist[i];
                maxDist[i] = m;
            }
            else if (m > maxDist2[i]) {
                maxDist2[i] = m;
            }
        },
        [&](node w, node z, count dist){
            if (sssp[i].distance(w) == sssp[i].distance(z)+1) { //update the shortest paths
                npaths[i][w] += npaths[i][z];
            }
        },
        [&](node w, count dist){ //node w's distance has become infinity
            DEBUG("distance of node ",w," has become infinity");
            vis[w]--;
            if (vis[w] == 0) { //if now no source is spanning w, we add it to U
                U.push(w);
            }
        });
        DEBUG("Done updating distance");
        if (maxDist[i] + maxDist2[i] > new_vd) {
            new_vd = maxDist[i] + maxDist2[i];
            INFO("VD has increased: before ", vd, " now ", new_vd);
        }
        if (sssp[i].modified()) {
            // subtract contributions to nodes in the old sampled path
            for (node z: sampledPaths[i]) {
                scoreData[z] -= 1 / (double) r;
            }
            if (sssp[i].distance(v[i]) < infDist) {
                // sample a new shortest path
                sampledPaths[i].clear();
                node t = v[i];
                while (t != u[i])  {
                    // sample z in P_u(t) with probability sigma_uz / sigma_us
                    std::vector<std::pair<node, double>> choices;
                    G.forEdgesOf(t, [&](node t, node z, edgeweight w){
                        if (Aux::NumericTools::logically_equal(sssp[i].distance(t), sssp[i].distance(z) + w)) {
                            // workaround for integer overflow in large graphs
                            double weight = npaths[i][z] / double(npaths[i][t]);
                            choices.emplace_back(z, weight);
                        }
                    });
                    assert (choices.size() > 0); // this should fail only if the graph is not connected
                    node z = Aux::Random::weightedChoice(choices);
                    assert (z <= G.upperNodeIdBound());
                    if (z != u[i]) {
                        scoreData[z] += 1 / (double) r;
                        sampledPaths[i].push_back(z);
                    }
                    t = z;
                }
            }

        }

    }
    INFO("Updating the additional r2 SSSPs");
    //now we update the additional r2 SSSPs
    for (node i = 0; i < r2; i++) {
        assert(sssp2.size() == r2);
        assert(maxDist.size() >= r + r2);
        assert(maxDist2.size() >= r + r2);
        sssp2[i].update(batch,
        [&](node w, count dist, count m){
            if(dist == infDist) { //node w has become reachable
                vis[w]++;
            }
                //calculate the maxDist
            if(m > maxDist[i+r]) {
                maxDist2[i+r] = maxDist[i+r];
                maxDist[i+r] = m;
            }
            else if (m > maxDist2[i+r]) {
                maxDist2[i+r] = m;
            }
        },
        [&](node w, node z, count dist){
        },
        [&](node w, count dist){ //node w's distance has become infinity
            vis[w]--;
            if (vis[w] == 0) //if now no source is spanning w, we add it to U
                U.push(w);
        });
        if (maxDist[i+r] + maxDist2[i+r] > new_vd) {
            new_vd = maxDist[i+r] + maxDist2[i+r];
        }
    }
    INFO("Sampling new additional SSSPs");
    while(!U.empty()) {
        node s = U.front();
        U.pop();
        if (vis[s] > 0) {
            continue;
        }
        count max1 = 0, max2 = 0;
        BFSvisit bfs(G,s);
        bfs.run([&](node w, count dist){
            vis[w]++;
            if(dist > max1) {
                max2 = max1;
                max1 = dist;
            }
            else if (dist > max2) {
                max2 = dist;
            }
        },
        [&](node w, node z, count dist_w, count dist_z, bool visited){
        });
        if (max1 + max2 > new_vd) {
            new_vd = max1 + max2;
        }
        sssp2.push_back(bfs);
        maxDist.push_back(max1);
        maxDist2.push_back(max2);
        r2 ++;
    }
    if (new_vd != vd) {
        assert(new_vd > vd);
        INFO("new vd: ", new_vd, ", old vd:", vd);
        vd = new_vd;
        count new_r = ceil((c / (epsilon * epsilon)) * (floor(log2(vd - 2)) + 1 + log(1 / delta)));
        if (new_r > r) {
            sampleNewPaths(r, new_r);
            //multiply all scores by r/new_r
            G.forNodes([&](node n) {
                scoreData[n] *= r/double(new_r);
            });
            r = new_r;
        }
    }
}

} /* namespace NetworKit */
