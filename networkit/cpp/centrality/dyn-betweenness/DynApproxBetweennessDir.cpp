/*
 * DynApproxBetweennessDir.cpp
 *
 *  Created on: 31.07.2014
 *      Author: ebergamini
 */

<<<<<<< Updated upstream
<<<<<<< HEAD
=======
// directed/undirected, unweighted

>>>>>>> Stashed changes
#include "DynApproxBetweennessDir.h"
#include "../auxiliary/Random.h"
#include "../properties/Diameter.h"
#include "../graph/Sampling.h"
#include "../graph/DynDijkstra.h"
#include "../graph/DynBFS.h"
#include "../auxiliary/Log.h"
#include "../properties/StronglyConnectedComponents.h"
#include "../auxiliary/NumericTools.h"
#include "../graph/BFSvisit.h"
#include <ctime>
#include <math.h>
=======
#include <networkit/centrality/DynApproxBetweennessDir.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/distance/Diameter.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/distance/DynDijkstra.hpp>
#include <networkit/distance/DynBFS.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>

#include <ctime>
#include <cmath>
>>>>>>> 7ca6f73b0... Fix includes and first steps at compiling files (not compiling yet)

namespace NetworKit {

DynApproxBetweennessDir::DynApproxBetweennessDir(const Graph& G, double epsilon, double delta, bool storePredecessors) : Centrality(G, true),
storePreds(storePredecessors), epsilon(epsilon), delta(delta) {

}


count DynApproxBetweennessDir::getNumberOfSamples() {
    return r;
}


inline count DynApproxBetweennessDir::DFS(node u, count j, Graph & sccDAG, std::vector<bool>& marked, std::vector<node>& sorted) {
    marked[u] = true;
    sccDAG.forNeighborsOf(u, [&](node v) {
        if (!marked[v]) {
            marked[v] = true;
            j = DFS(v, j, sccDAG, marked, sorted);
        }
    });
    j --;
    sorted[j] = u;
    return j;
}

count DynApproxBetweennessDir::computeVDdirected() {
  clock_t begin = clock();
  StronglyConnectedComponents scc(G);
  scc.run();
  clock_t end = clock();
  INFO("Time to run dyn conn comp: ", end-begin);
  begin = clock();
  count sccNum = scc.numberOfComponents();
  INFO("Number of nodes: ", G.upperNodeIdBound(), ", number of components: ", sccNum);


  Graph sccDAG(sccNum, true, true);
  //Partition sccPart = scc.getPartition();
  // std::map<index, count> sccSizes = sccPart.subsetSizeMap();
  count maxSize = 0;
  G.forEdges([&](node u, node v){
    if (scc.componentOfNode(u)!=scc.componentOfNode(v)) {
      if(!sccDAG.hasEdge(scc.componentOfNode(u)-1,scc.componentOfNode(v)-1)) {
        sccDAG.addEdge(scc.componentOfNode(u)-1,scc.componentOfNode(v)-1); // maybe sccSizes[scc.componentOfNode(v)] ?
      }
      // if (sccSizes[scc.componentOfNode(v)] > maxSize) {
      //   maxSize = sccSizes[scc.componentOfNode(v)];
      // }
    }
  });
  end = clock();
  INFO("Time to create scc DAG: ", end-begin);
  // INFO("Max component size: ", maxSize);
  // compute topological sorting of the connected components DAG
  begin = clock();
  node r; // root of the DAG
  count j = sccNum;
  std::vector<node> sorted(sccNum);
  std::vector<bool> marked(sccNum);


  sccDAG.forNodes([&](node c){
      if (!marked[c] && sccDAG.degreeIn(c) == 0) {
        marked[c] = true;
        j = DFS(c, j, sccDAG, marked, sorted);
      }
    });
 //  sccDAG.forEdges([&](node u, node v, edgeweight w){
 //    INFO("Edge ",u,",",v,", weight ",w);
 //  });
 //  for (count i = 0; i < sccNum; i ++) {
 //   	INFO("Comp [",i,"] = ", sorted[i]);
 // }

 end = clock();
 INFO("Time for topological sorting: ", end-begin);
 begin = clock();
 std::vector<bool> marked1(sccNum);
 marked1.resize(G.upperNodeIdBound(), false);
 std::vector<bool> marked2(sccNum);
 marked2.resize(G.upperNodeIdBound(), false);
 std::vector<node> sources(sccNum, 0);
 G.forNodes([&](node u){
   count x = scc.componentOfNode(u)-1;
   if (sources[x] == 0) {
     sources[x] = u;
   }
 });
  // we compute the bound on VD for each SCC
  std::vector<count> vdapproxs(sccNum);
  std::queue<node> q, qNext;
  for (count i = 0 ; i< sccNum; i++) {
    //count source = * (sccPart.getMembers(sorted[i]+1).begin());
    count source = sources[sorted[i]];
    //count source = *members.begin(); // we get a node from the set
    // we run one forward and one backward BFS from source, only within the SCC of source
    // forward BFS
    marked1[source] = true;
    count dist1 = 0;
    q.push(source);
    do {
        node u = q.front();
        q.pop();
        // apply function
        G.forNeighborsOf(u, [&](node v) {
            if (!marked1[v] && scc.componentOfNode(u)==scc.componentOfNode(v)) {
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
        if (!marked2[v] && scc.componentOfNode(u)==scc.componentOfNode(v)) {
          qNext.push(v);
          marked2[v] = true;
        }
      });
      if (q.empty() && !qNext.empty()) {
        q.swap(qNext);
        ++dist2;
      }
    } while (!q.empty());
    DEBUG("Component ", sorted[i]+1, ". Dist1 = ", dist1, ", dist2 = ", dist2, ", computed form node ", source);
    vdapproxs[sorted[i]] = dist1 + dist2 + 1; // source node is not taken into account in any of the two approximations
    if (vdapproxs[sorted[i]] <= 0) { // for singletons
      vdapproxs[sorted[i]] = 1;
    }
    // if (sccSizes[sorted[i]+1] < vdapproxs[sorted[i]]) { // can this actually ever happen?
    //   vdapproxs[sorted[i]] = sccSizes[sorted[i]+1];
    // }
  }
  end = clock();
  INFO("Time for VD approximations: ", end-begin);
  // starting from the bottom, we compute all the cumulative vd approximations
  begin = clock();
  count vd = vdapproxs[0];
  std::vector<count> vdSuccessors(sccNum, 0);
  count vd_comp;
  for (int i = sccNum-1; i >=0; i--) {
  //  INFO("Extracted component ", sorted[i]+1, ", i = ", i);
    node c = sorted[i]; // should be sorted[i]+1
    if (sccDAG.degreeOut(c) > 0) {
      vdapproxs[c] += vdSuccessors[c];
    }
    if (vdapproxs[c] > vd) {
      vd = vdapproxs[c];
      vd_comp = c+1;
    }
    sccDAG.forInEdgesOf(c, [&](node c, node c_pred, edgeweight w){
      if (vdSuccessors[c_pred] < vdapproxs[c])
        vdSuccessors[c_pred] = vdapproxs[c];
    });
  }
  end = clock();
  INFO("Time for finding maximum: ", end-begin);
  // INFO("size of the initial component: ", sccSizes[vd_comp]);
  // std::set<index> members = sccPart.getMembers(vd_comp);
  // std::set<index>::iterator it;
  // count maxdist = 0;
  // for (it = members.begin(); it != members.end(); ++it){
  //   count source = *it;
  //   G.BFSfrom(source, [&](node v, count dist) {
  //     if (dist > maxdist) {
  //       maxdist = dist;
  //     }
  //   });
  // }
  //
  // INFO("MAXDIST = ", maxdist);

  // for (count i = 0; i < sccNum; i ++) {
  //   INFO("vdapproxs [",i,"] = ", vdapproxs[i]);
  // }
  INFO("VD = ", vd);
  return vd;
}

void DynApproxBetweennessDir::sampleNewPaths(count start, count end) {
    for (count i = start; i < end; i++) {
        DEBUG("sample ", i);
        // sample random node pair
        node u1, u2;
<<<<<<< HEAD
        u1 = Sampling::randomNode(G);
        do {
            u2 = Sampling::randomNode(G);
=======
        u1 = GraphTools::randomNode(G);
        do {
            u2 = GraphTools::randomNode(G);
>>>>>>> 7ca6f73b0... Fix includes and first steps at compiling files (not compiling yet)
        } while (u1 == u2);
        u.push_back(u1);
        v.push_back(u2);
        if (G.isWeighted()) {
            return; //TODO weighted graphs!
        //    sssp[i].reset(new DynDijkstra(G, u[i], storePreds));
        } else {
<<<<<<< HEAD
            BFSvisit bfs(G, u[i]);
=======
			DynBFS bfs(G, u[i], storePreds);
>>>>>>> 7ca6f73b0... Fix includes and first steps at compiling files (not compiling yet)
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
                G.forInEdgesOf(t, [&](node t, node z, edgeweight w){
                    if (Aux::NumericTools::logically_equal(sssp[i].distance(t), sssp[i].distance(z) + w)) {
                        // workaround for integer overflow in large graphs
                        double weight = npaths[i][z] / double(npaths[i][t]);
                        assert(weight != 0);
                        choices.emplace_back(z, weight);
                    }

                });
<<<<<<< HEAD
                assert (choices.size() > 0);
=======
                assert (!choices.empty() > 0);
>>>>>>> 7ca6f73b0... Fix includes and first steps at compiling files (not compiling yet)
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

void DynApproxBetweennessDir::run() {
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

  count vd;
  if (G.isDirected())
  {
    vd = computeVDdirected();
  } else {
<<<<<<< HEAD
    vd = Diameter::estimatedVertexDiameterPedantic(G);
=======
    Diameter diam(G, DiameterAlgo::estimatedPedantic);
    diam.run();
    vd = diam.getDiameter().first;
>>>>>>> 7ca6f73b0... Fix includes and first steps at compiling files (not compiling yet)
  }

  INFO("estimated diameter: ", vd);
  r = ceil((c / (epsilon * epsilon)) * (floor(log2(vd - 2)) + 1 + log(1 / delta)));
  INFO("taking ", r, " path samples");
  sssp.clear();

  sampleNewPaths(0, r);
  hasRun = true;

}

void DynApproxBetweennessDir::update(const std::vector<GraphEvent>& batch) {
  std::queue<node> U;
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
                  DEBUG("npaths = ", npaths[i][t]);
                  DEBUG("distance = ", sssp[i].distance(t));
                  DEBUG("Node ", t, ", source = ", u[i], ". Is there edge (s,t)? ", G.hasEdge(u[i], t));
                  G.forInEdgesOf(t, [&](node t, node z, edgeweight w){
                      if (Aux::NumericTools::logically_equal(sssp[i].distance(t), sssp[i].distance(z) + w)) {
                          // workaround for integer overflow in large graphs
                          double weight = npaths[i][z] / double(npaths[i][t]);
                          choices.emplace_back(z, weight);
                      }
                  });
<<<<<<< HEAD
                  assert (choices.size() > 0); // this should fail only if the graph is not connected
=======
                  assert (!choices.empty()); // this should fail only if the graph is not connected
>>>>>>> 7ca6f73b0... Fix includes and first steps at compiling files (not compiling yet)
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
  count new_vd = computeVDdirected();
  count new_r = ceil((c / (epsilon * epsilon)) * (floor(log2(new_vd - 2)) + 1 + log(1 / delta)));
  if (new_r > r) {
      sampleNewPaths(r, new_r);
      //multiply all scores by r/new_r
      G.forNodes([&](node n) {
          scoreData[n] *= r/double(new_r);
      });
      r = new_r;
  }
}

} /* namespace NetworKit */
