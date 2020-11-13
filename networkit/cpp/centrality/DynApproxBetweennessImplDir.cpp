/*
 * DynApproxBetweennssImplDir.cpp
 *
 * Author: Charmaine Ndolo <Charmaine.Ndolo@hu-berlin.de>
 */

// networkit-format
#include "DynApproxBetweennessImplDir.hpp"
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/centrality/Centrality.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/distance/DynBFS.hpp>
#include <networkit/distance/DynDijkstra.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

DynApproxBetweennessImplDir::DynApproxBetweennessImplDir(const Graph &G, double epsilon,
                                                         double delta, bool storePredecessors)
    : Centrality(G, true), epsilon(epsilon), delta(delta), storePreds(storePredecessors) {}

inline count DynApproxBetweennessImplDir::DFS(node u, count j, Graph &sccDAG,
                                              std::vector<bool> &marked,
                                              std::vector<node> &sorted) {
    marked[u] = true;
    sccDAG.forNeighborsOf(u, [&](node v) {
        if (!marked[v]) {
            marked[v] = true;
            j = DFS(v, j, sccDAG, marked, sorted);
        }
    });
    j--;
    sorted[j] = u;
    return j;
}

count DynApproxBetweennessImplDir::computeVDdirected() {
    clock_t begin = clock();
    StronglyConnectedComponents scc(G);
    scc.run();
    clock_t end = clock();
    INFO("Time to run dyn conn comp: ", end - begin);
    begin = clock();
    count sccNum = scc.numberOfComponents();
    INFO("Number of nodes: ", G.upperNodeIdBound(), ", number of components: ", sccNum);
    Graph sccDAG(sccNum, true, true);
    count maxSize = 0;
    G.forEdges([&](node u, node v) {
        if (scc.componentOfNode(u) != scc.componentOfNode(v)) {
            if (!sccDAG.hasEdge(scc.componentOfNode(u) - 1, scc.componentOfNode(v) - 1)) {
                sccDAG.addEdge(scc.componentOfNode(u) - 1,
                               scc.componentOfNode(v)
                                   - 1); // maybe sccSizes[scc.componentOfNode(v)] ?
            }
        }
    });
    end = clock();
    INFO("Time to create scc DAG: ", end - begin);
    begin = clock();
    node r; // root of the DAG
    count j = sccNum;
    std::vector<node> sorted(sccNum);
    std::vector<bool> marked(sccNum);

    sccDAG.forNodes([&](node c) {
        if (!marked[c] && sccDAG.degreeIn(c) == 0) {
            marked[c] = true;
            j = DFS(c, j, sccDAG, marked, sorted);
        }
    });

    end = clock();
    INFO("Time for topological sorting: ", end - begin);
    begin = clock();
    std::vector<bool> marked1(sccNum);
    marked1.resize(G.upperNodeIdBound(), false);
    std::vector<bool> marked2(sccNum);
    marked2.resize(G.upperNodeIdBound(), false);
    std::vector<node> sources(sccNum, 0);
    G.forNodes([&](node u) {
        count x = scc.componentOfNode(u) - 1;
        if (sources[x] == 0) {
            sources[x] = u;
        }
    });
    // we compute the bound on VD for each SCC
    std::vector<count> vdapproxs(sccNum);
    std::queue<node> q, qNext;
    for (count i = 0; i < sccNum; i++) {
        count source = sources[sorted[i]];
        // forward BFS
        marked1[source] = true;
        count dist1 = 0;
        q.push(source);
        do {
            node u = q.front();
            q.pop();
            // apply function
            G.forNeighborsOf(u, [&](node v) {
                if (!marked1[v] && scc.componentOfNode(u) == scc.componentOfNode(v)) {
                    qNext.push(v);
                    marked1[v] = true;
                }
            });
            if (q.empty() && !qNext.empty()) {
                q.swap(qNext);
                ++dist1;
            }
        } while (!q.empty());

        // backward BFS
        marked2[source] = true;
        count dist2 = 0;
        q.push(source);
        do {
            node u = q.front();
            q.pop();
            // apply function
            G.forInNeighborsOf(u, [&](node v) {
                if (!marked2[v] && scc.componentOfNode(u) == scc.componentOfNode(v)) {
                    qNext.push(v);
                    marked2[v] = true;
                }
            });
            if (q.empty() && !qNext.empty()) {
                q.swap(qNext);
                ++dist2;
            }
        } while (!q.empty());
        DEBUG("Component ", sorted[i] + 1, ". Dist1 = ", dist1, ", dist2 = ", dist2,
              ", computed form node ", source);
        vdapproxs[sorted[i]] =
            dist1 + dist2
            + 1; // source node is not taken into account in any of the two approximations
        if (vdapproxs[sorted[i]] <= 0) { // for singletons
            vdapproxs[sorted[i]] = 1;
        }
    }
    end = clock();
    INFO("Time for VD approximations: ", end - begin);
    // starting from the bottom, we compute all the cumulative vd approximations
    begin = clock();
    count vd = vdapproxs[0];
    std::vector<count> vdSuccessors(sccNum, 0);
    count vd_comp;
    for (int i = sccNum - 1; i >= 0; i--) {
        node c = sorted[i]; // should be sorted[i]+1
        if (sccDAG.degreeOut(c) > 0) {
            vdapproxs[c] += vdSuccessors[c];
        }
        if (vdapproxs[c] > vd) {
            vd = vdapproxs[c];
            vd_comp = c + 1;
        }
        sccDAG.forInEdgesOf(c, [&](node c, node c_pred, edgeweight w) {
            if (vdSuccessors[c_pred] < vdapproxs[c])
                vdSuccessors[c_pred] = vdapproxs[c];
        });
    }
    end = clock();
    INFO("Time for finding maximum: ", end - begin);
    INFO("VD = ", vd);
    return vd;
}

count DynApproxBetweennessImplDir::getNumberOfSamples() {
    return r;
}

void DynApproxBetweennessImplDir::run() {
    scoreData.clear();
    scoreData.resize(G.upperNodeIdBound());
    u.clear();
    v.clear();
    sampledPaths.clear();
    vis.clear();
    vis.resize(G.upperNodeIdBound());
    npaths.clear();
    // maxDist.clear();
    // maxDist2.clear();

    count vd;
    if (G.isDirected()) {
        vd = computeVDdirected();
    } else {
        Diameter diam(G, DiameterAlgo::estimatedPedantic);
        diam.run();
        vd = diam.getDiameter().first;
    }

    INFO("estimated diameter: ", vd);
    r = ceil((c / (epsilon * epsilon)) * (floor(log2(vd - 2)) + 1 + log(1 / delta)));
    INFO("taking ", r, " path samples");
    sssp.clear();

    sampleNewPaths(0, r);
    hasRun = true;
}

void DynApproxBetweennessImplDir::update(GraphEvent e) {
    std::vector<GraphEvent> batch(1);
    batch[0] = e;
    updateBatch(batch);
}

void DynApproxBetweennessImplDir::updateBatch(const std::vector<GraphEvent> &batch) {
    std::queue<node> U;
    INFO("Updating the r paths");
    assert(sssp.size() == r);
    for (node i = 0; i < r; i++) {
        DEBUG("Updating distance of source node number ", i, " out of ", r);
        assert(npaths.size() == r);
        assert(vis.size() == G.upperNodeIdBound());
        assert(maxDist.size() >= r);
        assert(maxDist2.size() >= r);
        sssp[i]->updateBatch(
            batch,
            [&](node w, count dist, count m) {
                DEBUG("Updating node ", w, ", source node: ", u[i], ", distance: ", dist,
                      ", m: ", m);
                assert(npaths[i].size() == G.upperNodeIdBound());
                npaths[i][w] = 0;
                if (dist == infDist) { // node w has become reachable
                    vis[w]++;
                }
                // calculate the maxDist
                if (m > maxDist[i]) {
                    maxDist2[i] = maxDist[i];
                    maxDist[i] = m;
                } else if (m > maxDist2[i]) {
                    maxDist2[i] = m;
                }
            },
            [&](node w, node z, count dist) {
                if (sssp[i]->distance(w) == sssp[i]->distance(z) + 1) { // update the shortest paths
                    npaths[i][w] += npaths[i][z];
                }
            },
            [&](node w, count dist) { // node w's distance has become infinity
                DEBUG("distance of node ", w, " has become infinity");
                vis[w]--;
                if (vis[w] == 0) { // if now no source is spanning w, we add it to U
                    U.push(w);
                }
            });
        DEBUG("Done updating distance");
        if (sssp[i]->modified()) {
            // subtract contributions to nodes in the old sampled path
            for (node z : sampledPaths[i]) {
                scoreData[z] -= 1 / (double)r;
            }
            if (sssp[i]->distance(v[i]) < infDist) {
                // sample a new shortest path
                sampledPaths[i].clear();
                node t = v[i];
                while (t != u[i]) {
                    // sample z in P_u(t) with probability sigma_uz / sigma_us
                    std::vector<std::pair<node, double>> choices;
                    DEBUG("npaths = ", npaths[i][t]);
                    DEBUG("distance = ", sssp[i]->distance(t));
                    DEBUG("Node ", t, ", source = ", u[i], ". Is there edge (s,t)? ",
                          G.hasEdge(u[i], t));
                    G.forInEdgesOf(t, [&](node t, node z, edgeweight w) {
                        if (Aux::NumericTools::logically_equal(sssp[i]->distance(t),
                                                               sssp[i]->distance(z) + w)) {
                            // workaround for integer overflow in large graphs
                            double weight = npaths[i][z] / double(npaths[i][t]);
                            choices.emplace_back(z, weight);
                        }
                    });
                    assert(choices.size()
                           > 0); // this should fail only if the graph is not connected
                    node z = Aux::Random::weightedChoice(choices);
                    assert(z <= G.upperNodeIdBound());
                    if (z != u[i]) {
                        scoreData[z] += 1 / (double)r;
                        sampledPaths[i].push_back(z);
                    }
                    t = z;
                }
            }
        }
    }
    count new_vd = computeVDdirected();
    count new_r = ceil((c / (epsilon * epsilon)) * (floor(log2(new_vd - 2)) + 1 + log(1 / delta)));
    if (new_r > r) {
        sampleNewPaths(r, new_r);
        // multiply all scores by r/new_r
        G.forNodes([&](node n) { scoreData[n] *= r / double(new_r); });
        r = new_r;
    }
}

void DynApproxBetweennessImplDir::sampleNewPaths(count start, count end) {
    for (count i = start; i < end; i++) {
        DEBUG("sample ", i);
        // sample random node pair
        node u1, u2;
        u1 = GraphTools::randomNode(G);
        do {
            u2 = GraphTools::randomNode(G);
        } while (u1 == u2);
        u.push_back(u1);
        v.push_back(u2);
        if (G.isWeighted()) {
            return;
        } else {
            sssp[i].reset(new DynBFS(G, u[i], storePreds));
        }
        DEBUG("running shortest path algorithm for node ", u[i]);
        std::vector<count> sigma(G.upperNodeIdBound(), 0);
        sigma[u[i]] = 1;
        count max1 = 0, max2 = 0;
        sssp[i]->run();
        npaths.push_back(sigma);
        // maxDist.push_back(max1);
        // maxDist2.push_back(max2);
        std::vector<node> path;
        if (sssp[i]->distance(v[i]) < infDist) { // at least one path between {u, v} exists
            DEBUG("updating estimate for path ", u[i], " <-> ", v[i]);
            // random path sampling and estimation update
            sampledPaths[i].clear();
            node t = v[i];
            while (t != u[i]) {
                // sample z in P_u(t) with probability sigma_uz / sigma_us
                std::vector<std::pair<node, double>> choices;
                if (storePreds) {
                    for (node z : sssp[i]->previous[t]) {
                        // workaround for integer overflow in large graphs
                        bigfloat tmp = sssp[i]->numberOfPaths(z) / sssp[i]->numberOfPaths(t);
                        double weight;
                        tmp.ToDouble(weight);

                        choices.emplace_back(z, weight); // sigma_uz / sigma_us
                    }
                } else {
                    G.forInEdgesOf(t, [&](node t, node z, edgeweight w) {
                        if (Aux::NumericTools::logically_equal(sssp[i]->distances[t],
                                                               sssp[i]->distances[z] + w)) {
                            // workaround for integer overflow in large graphs
                            bigfloat tmp = sssp[i]->numberOfPaths(z) / sssp[i]->numberOfPaths(t);
                            double weight;
                            tmp.ToDouble(weight);

                            choices.emplace_back(z, weight);
                        }
                    });
                }
                DEBUG("Node: ", t);
                DEBUG("Source: ", u[i]);
                assert(!choices.empty());
                node z = Aux::Random::weightedChoice(choices);
                assert(z <= G.upperNodeIdBound());
                if (z != u[i]) {
                    scoreData[z] += 1 / (double)r;
                    sampledPaths[i].push_back(z);
                }
                t = z;
            }
        }
    }
}
} // namespace NetworKit
