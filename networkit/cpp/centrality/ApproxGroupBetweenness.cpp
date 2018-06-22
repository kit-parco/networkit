/*
 * ApproxGroupBetweenness.h
 *
 *  Created on: 13.03.2018
 *      Author: Marvin Pogoda
 */

#include "../auxiliary/BucketPQ.h"
#include "../auxiliary/Random.h"
#include "../distance/BFS.h"
#include "../distance/SSSP.h"
#include "ApproxGroupBetweenness.h"

#include <math.h>
#include <omp.h>

namespace NetworKit {

ApproxGroupBetweenness::ApproxGroupBetweenness(const Graph &G,
                                               const count groupSize,
                                               const double epsilon)
    : G(G), groupSize(groupSize), epsilon(epsilon) {
	hasRun = false;
}

void ApproxGroupBetweenness::run() {
	// Create data structures for the hypergraph.
	std::vector<int64_t> bucketInitializer(G.numberOfNodes());
	std::vector<std::vector<node>> incidencyList(G.numberOfNodes());
	std::vector<count> hyperEdges;
	count samples = groupSize * log(G.numberOfNodes()) / (epsilon * epsilon);
	std::vector<std::vector<count>> hyperEdgesPerSample(samples);
	Aux::BucketPQ nodeDegrees(bucketInitializer, -samples, 1);

#pragma omp parallel for
	for (count l = 0; l < samples; ++l) {
		node s = G.randomNode();
		node t;
		do {
			t = G.randomNode();
		} while (s == t);

		BFS bfs(G, s, true, true, t);
		bfs.run();
		std::set<std::vector<node>> shortestPaths = bfs.getPaths(t);

		// If the selected nodes are in different connected components, the
		// hyperedge is an empty set. Chooseing nodes in different connected
		// components wont affect the algorithm. (See Mahmoody "Scalable Betweenness
		// Centrality Maximization via Sampling",page 4,Lemma 3,2016)
		if (shortestPaths.size() == 0) {
			continue;
		}

		// Uniformly select a shortest path
		std::set<std::vector<node>>::const_iterator iterator(shortestPaths.begin());
		advance(iterator, Aux::Random::integer(shortestPaths.size() - 1));
		std::vector<node> newHyperEdge = *iterator;

		// Insert the HyperEdge into hyperEdgePerSample
		for (auto n : newHyperEdge) {
			hyperEdgesPerSample[l].push_back(n);
		}
	}

	// Transfer edges from hyperEdgesPerSample to hyperEdges and prepare building
	// nodeDegrees.
	std::vector<count> tempDegrees(G.numberOfNodes());
	for (auto edge : hyperEdgesPerSample) {
		node hyperEdgeStart = hyperEdges.size();
		hyperEdges.push_back(edge.size());
		for (auto const &n : edge) {
			// Safe degrees temporary to minimize the use of the changeKey-option
			tempDegrees[n] -= 1;
			// Insert new Hyperedge into the hyperedge-list
			hyperEdges.push_back(n);
			// Update incidencecyList
			incidencyList[n].push_back(hyperEdgeStart);
		}
	}
	// Build nodeDegrees
	for (node i = 0; i < G.numberOfNodes(); i++)
		nodeDegrees.changeKey(tempDegrees[i], i);

	// Extract nodes with highest degrees.
	std::vector<count> degreeDecrease(G.numberOfNodes());
	node v;
	count degree;
	for (count j = 0; j < groupSize; j++) {
		// Lazy-Queue-Update tries to reduce the number of updates,that need to be
		// done to get the next maximum degree. We only need to update the
		// max-degree-node as best case. The worst case is updateing every node.
		std::pair<count, node> elem = nodeDegrees.extractMin();
		v = elem.second;
		degree = elem.first;
		// degreeDecreased[v]=0 -> v has max degree
		while (degreeDecrease[v] > 0) {
			// Update.
			degree += degreeDecrease[v];
			degreeDecrease[v] = 0;
			// Check if v still has max degree
			nodeDegrees.changeKey(degree, v);
			elem = nodeDegrees.extractMin();
			v = elem.second;
			degree = elem.first;
		}

		// Upgrade temporal degreeDecrease to minimize number of accesses of
		// BucketPQ
		for (auto hyperEdge : incidencyList[v]) {
			count start = hyperEdge + 1;
			count end = start + hyperEdges[hyperEdge];
			for (count i = start; i < end; i++) {
				if (hyperEdges[i] != v) {
					degreeDecrease[hyperEdges[i]] += 1;
				}
			}
			//"Deleting" Hyperedge
			hyperEdges[start - 1] = 0;
		}
		maxGroup.push_back(v);
	}

	hasRun = true;
}

} /* namespace NetworKit */