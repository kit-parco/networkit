/*
 * DynApproxBetweenness2W.cpp
 *
 *  Created on: 05.03.2015
 *      Author: ebergamini
 */

#include "DynApproxBetweenness2W.h"
#include "../auxiliary/Random.h"
#include "../properties/Diameter.h"
#include "../graph/Sampling.h"
#include "../graph/DynDijkstra.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/NumericTools.h"
#include <queue>
#include <math.h>


namespace NetworKit {

DynApproxBetweenness2W::DynApproxBetweenness2W(const Graph& G, double epsilon, double delta, bool storePredecessors) : Centrality(G, true),
storePreds(storePredecessors), epsilon(epsilon), delta(delta) {

}


count DynApproxBetweenness2W::getNumberOfSamples() {
    return r;
}


void DynApproxBetweenness2W::run() {
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
    wmin.clear();
    compSize.clear();

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
        DijkstraVisit dij(G,s);
        count max1 = 0, max2 = 0;
        count thisCompSize = 0;
        edgeweight wminLocal = infDist;
        dij.run([&](node w, count dist){
            vis[w]++;
            thisCompSize ++;
            if (dist > max1) {
                max2 = max1;
                max1 = dist;
            } else if (dist > max2) {
                max2 = dist;
            }
        },
        [&](node w, node z, count dist_w, count dist_z, edgeweight ew){
            if (ew < wminLocal) {
                wminLocal = ew;
            }
        });
        sssp2.push_back(dij);
        maxDist.push_back(max1);
        maxDist2.push_back(max2);
        wmin.push_back(wminLocal);
        compSize.push_back(thisCompSize);
        r2 ++;
    }
    INFO("SampledPaths size: ", sampledPaths.size());
    assert(sampledPaths.size() == r);
    assert(npaths.size() == r);
    assert(maxDist.size() == r+r2);
    assert(maxDist2.size() == r+r2);
    hasRun = true;
}


void DynApproxBetweenness2W::sampleNewPaths(count start, count end) {
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
        DijkstraVisit dij(G, u[i]);
        sssp.push_back(dij);
        DEBUG("running shortest path algorithm for node ", u[i]);
        std::vector<count> sigma(G.upperNodeIdBound(), 0);
        sigma[u[i]] = 1;
        count max1 = 0, max2 = 0;
        count thisCompSize = 0;
        edgeweight wminLocal = infDist;
        assert(sssp.size() == i+1);
        sssp[i].run([&](node w, count dist) {
            if (dist < infDist) { //TODO: this is a difference between weighted and unweighted!
                vis[w]++;
                thisCompSize ++;
                if (dist > max1) {
                    max2 = max1;
                    max1 = dist;
                } else if (dist > max2) {
                    max2 = dist;
                }
            }
        },
        //we build sigma
        [&](node w, node z, edgeweight dist_w, edgeweight dist_z, edgeweight ew) {
            if (ew < wminLocal) {
                wminLocal = ew;
            }
            if (dist_z > dist_w + ew) { //TODO this is different from the unweighted case
                sigma[z] = sigma[w];
            } else if (dist_z == dist_w + ew) {
                sigma[z] += sigma[w];
            }
            assert(sigma[z] > 0 || dist_z == infDist);
        });
        npaths.push_back(sigma);
        maxDist.push_back(max1);
        maxDist2.push_back(max2);
        wmin.push_back(wminLocal);
        compSize.push_back(thisCompSize);
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

void DynApproxBetweenness2W::update(const std::vector<GraphEvent>& batch) {
    std::queue<node> U;
    count new_vd = vd;
    DEBUG("Updating the r paths");
    for (node i = 0; i < r; i++) {
        DEBUG("Updating distance of source node number ",i," out of ",r);
        sssp[i].update(batch,
        [&](node w, count dist, count m){
            npaths[i][w] = 0;
            if(dist == infDist) { //node w has become reachable
                vis[w]++;
                compSize[i]++;
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
            if (sssp[i].distance(w) == sssp[i].distance(z)+G.weight(w,z)) { //update the shortest paths
                npaths[i][w] += npaths[i][z];
            }
            if (G.weight(w,z) < wmin[i]) {
                wmin[i] = G.weight(w,z);
            }

        },
        [&](node w, count dist){ //node w's distance has become infinity
            vis[w]--;
            compSize[i]--;
            if (vis[w] == 0) { //if now no source is spanning w, we add it to U
                U.push(w);
            }
        });
        DEBUG("Done updating distance");
        edgeweight localVD;
        if (compSize[i] < 1 + (maxDist[i] + maxDist2[i])/wmin[i])
            localVD = compSize[i];
        else
            localVD = 1 + (maxDist[i] + maxDist2[i])/wmin[i];
        if (localVD > new_vd) {
            new_vd = localVD;
            INFO("VD has increased: before ", vd, " now ", new_vd);
        }
        //updating the path
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
                    assert (choices.size() > 0);
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
        sssp2[i].update(batch,
        [&](node w, count dist, count m){
            if(dist == infDist) { //node w has become reachable
                compSize[i+r]++;
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
            if (G.weight(w,z) < wmin[i+r]) {
                wmin[i+r] = G.weight(w,z);
            }
        },
        [&](node w, count dist){ //node w's distance has become infinity
            vis[w]--;
            compSize[i+r]--;
            if (vis[w] == 0) //if now no source is spanning w, we add it to U
                U.push(w);
        });
        edgeweight localVD;
        if (compSize[i] < 1 + (maxDist[i+r] + maxDist2[i+r])/wmin[i+r])
            localVD = compSize[i+r];
        else
            localVD = 1 + (maxDist[i+r] + maxDist2[i+r])/wmin[i+r];
        if (localVD > new_vd) {
            new_vd = localVD;
            INFO("VD has increased: before ", vd, " now ", new_vd);
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
        edgeweight wminLocal = infDist;
        count compSizeLocal = 0;
        DijkstraVisit dij(G,s);
        dij.run([&](node w, count dist){
            vis[w]++;
            compSizeLocal ++;
            if(dist > max1) {
                max2 = max1;
                max1 = dist;
            }
            else if (dist > max2) {
                max2 = dist;
            }
        },
        [&](node w, node z, count dist_w, count dist_z, edgeweight ew){
            if (ew < wminLocal)
                wminLocal = ew;
        });
        sssp2.push_back(dij);
        maxDist.push_back(max1);
        maxDist2.push_back(max2);
        r2 ++;
        edgeweight localVD;
        if (compSizeLocal < 1 + (max1 + max2)/wminLocal)
            localVD = compSizeLocal;
        else
            localVD = 1 + (max1 + max2)/wminLocal;
        if (localVD > new_vd) {
            new_vd = localVD;
            INFO("VD has increased: before ", vd, " now ", new_vd);
        }
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
