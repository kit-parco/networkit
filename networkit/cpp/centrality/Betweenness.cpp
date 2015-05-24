/*
 * Betweenness.cpp
 *
 *  Created on: 29.07.2014
 *      Author: cls, ebergamini
 */

#include <stack>
#include <queue>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Betweenness.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Log.h"
#include "../graph/SSSP.h"
#include "../graph/Dijkstra.h"
#include "../graph/BFS.h"

namespace NetworKit {

Betweenness::Betweenness(const Graph& G, bool normalized, bool computeEdgeCentrality) : Centrality(G, normalized, computeEdgeCentrality) {

}

void Betweenness::run() {
	count z = G.upperNodeIdBound();
	scoreData.clear();
	scoreData.resize(z);
	if (computeEdgeCentrality) {
		count z2 = G.upperEdgeIdBound();
		edgeScoreData.clear();
		edgeScoreData.resize(z2);
	}

	// thread-local scores for efficient parallelism
	count maxThreads = omp_get_max_threads();
	std::vector<std::vector<double> > scorePerThread(maxThreads, std::vector<double>(G.upperNodeIdBound()));
	INFO("score per thread: ", scorePerThread.size());
	INFO("G.upperEdgeIdBound(): ", G.upperEdgeIdBound());
	std::vector<std::vector<double> > edgeScorePerThread(maxThreads, std::vector<double>(G.upperEdgeIdBound()));
	INFO("edge score per thread: ", edgeScorePerThread.size());

	auto computeDependencies = [&](node s) {

		std::vector<double> dependency(z, 0.0);

		// run SSSP algorithm and keep track of everything
		std::unique_ptr<SSSP> sssp;
		if (G.isWeighted()) {
			sssp.reset(new Dijkstra(G, s, true, true));
		} else {
			sssp.reset(new BFS(G, s, true, true));
		}

		sssp->run();

		// compute dependencies for nodes in order of decreasing distance from s
		std::stack<node> stack = sssp->getStack();
		while (!stack.empty()) {
			node t = stack.top();
			stack.pop();
			for (node p : sssp->getPredecessors(t)) {
				// workaround for integer overflow in large graphs
				bigfloat tmp = sssp->numberOfPaths(p) / sssp->numberOfPaths(t);
				double weight;
				tmp.ToDouble(weight);
				double c= weight * (1 + dependency[t]);
				dependency[p] += c;
				if (computeEdgeCentrality) {
					edgeScorePerThread[omp_get_thread_num()][G.edgeId(p,t)] += c;
				}


			}
			if (t != s) {
				scorePerThread[omp_get_thread_num()][t] += dependency[t];
			}
		}
	};

	G.balancedParallelForNodes(computeDependencies);

	INFO("adding thread-local scores");
	// add up all thread-local values
	for (auto local : scorePerThread) {
		G.parallelForNodes([&](node v){
			scoreData[v] += local[v];
		});
	}
	if (computeEdgeCentrality) {
		for (auto local : edgeScorePerThread) {
			for (count i = 0; i < local.size(); i++) {
				edgeScoreData[i] += local[i];
			}
		}
	}
	if (normalized) {
		// divide by the number of possible pairs
		count n = G.numberOfNodes();
		count pairs = (n-2) * (n-1);
		count edges =  n    * (n-1);
		G.forNodes([&](node u){
			scoreData[u] = scoreData[u] / pairs;
		});
		if (computeEdgeCentrality) {
			for (count edge = 0; edge < edgeScoreData.size(); edge++) {
				edgeScoreData[edge] =  edgeScoreData[edge] / edges;
			}
		}
	}

	ran = true;
}

} /* namespace NetworKit */
