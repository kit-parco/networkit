/*
 * GraphBuilder.cpp
 *
 *  Created on: 15.07.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <omp.h>
#include <stdexcept>

#include <networkit/graph/GraphBuilder.hpp>

namespace NetworKit {

template<typename X, typename Y>
static void swapTransform(std::vector<X> &vecX, std::vector<Y> &vecY) {
	std::vector<X> temp;
	temp.swap(vecX);
	// Copy from vecX to vecY.
	vecX.resize(vecY.size());
	std::copy(vecY.begin(), vecY.end(), vecX.begin());
	// Copy from temp (= old vecX) to vecY.
	vecY.clear();
	vecY.resize(temp.size());
	std::copy(temp.begin(), temp.end(), vecY.begin());
}

// More specific overload of the function above.
template<typename X>
static void swapTransform(std::vector<X> &vecX, std::vector<X> &vecY) {
	vecX.swap(vecY);
}

GraphBuilder::GraphBuilder(count n, bool weighted, bool directed)
    : n(n), selfloops(0), name(""), weighted(weighted), directed(directed),
      outEdges(n), outEdgeWeights(weighted ? n : 0), inEdges(directed ? n : 0),
      inEdgeWeights((directed && weighted) ? n : 0) {}

void GraphBuilder::reset(count n) {
	this->n = n;
	selfloops = 0;
	name = "";
	outEdges.assign(n, std::vector<storednode>{});
	outEdgeWeights.assign(isWeighted() ? n : 0, std::vector<edgeweight>{}),
	    inEdges.assign(isDirected() ? n : 0, std::vector<storednode>{}),
	    inEdgeWeights.assign((isDirected() && isWeighted()) ? n : 0,
	                         std::vector<edgeweight>{});
}

index GraphBuilder::indexInOutEdgeArray(node u, node v) const {
	for (index i = 0; i < outEdges[u].size(); i++) {
		node x = outEdges[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

index GraphBuilder::indexInInEdgeArray(node u, node v) const {
	assert(isDirected());
	for (index i = 0; i < inEdges[u].size(); i++) {
		node x = inEdges[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

node GraphBuilder::addNode() {
	outEdges.push_back(std::vector<storednode>{});
	if (weighted) {
		outEdgeWeights.push_back(std::vector<edgeweight>{});
	}
	if (directed) {
		inEdges.push_back(std::vector<storednode>{});
		if (weighted) {
			inEdgeWeights.push_back(std::vector<edgeweight>{});
		}
	}
	return n++;
}

void GraphBuilder::addHalfOutEdge(node u, node v, edgeweight ew) {
	assert(indexInOutEdgeArray(u, v) == none);
	outEdges[u].push_back(v);
	if (weighted) {
		outEdgeWeights[u].push_back(ew);
	}
	if (u == v) {
#pragma omp atomic
		selfloops++;
	}
}

void GraphBuilder::addHalfInEdge(node u, node v, edgeweight ew) {
	assert(indexInInEdgeArray(u, v) == none);
	inEdges[u].push_back(v);
	if (weighted) {
		inEdgeWeights[u].push_back(ew);
	}
	if (u == v) {
#pragma omp atomic
		selfloops++;
	}
}

void GraphBuilder::swapNeighborhood(node u, std::vector<node> &neighbours,
                                    std::vector<edgeweight> &weights,
                                    bool selfloop) {
	if (weighted)
		assert(neighbours.size() == weights.size());
	swapTransform(outEdges[u], neighbours);
	if (weighted) {
		outEdgeWeights[u].swap(weights);
	}

	if (selfloop) {
#pragma omp atomic
		selfloops++;
	}
}
void GraphBuilder::setOutWeight(node u, node v, edgeweight ew) {
	assert(isWeighted());
	index vi = indexInOutEdgeArray(u, v);
	if (vi != none) {
		outEdgeWeights[u][vi] = ew;
	} else {
		addHalfOutEdge(u, v, ew);
	}
}

void GraphBuilder::setInWeight(node u, node v, edgeweight ew) {
	assert(isWeighted());
	assert(isDirected());
	index vi = indexInInEdgeArray(u, v);
	if (vi != none) {
		inEdgeWeights[u][vi] = ew;
	} else {
		addHalfInEdge(u, v, ew);
	}
}

void GraphBuilder::increaseOutWeight(node u, node v, edgeweight ew) {
	assert(isWeighted());
	index vi = indexInOutEdgeArray(u, v);
	if (vi != none) {
		outEdgeWeights[u][vi] += ew;
	} else {
		addHalfOutEdge(u, v, ew);
	}
}

void GraphBuilder::increaseInWeight(node u, node v, edgeweight ew) {
	assert(isWeighted());
	assert(isDirected());
	index vi = indexInInEdgeArray(u, v);
	if (vi != none) {
		inEdgeWeights[u][vi] += ew;
	} else {
		addHalfInEdge(u, v, ew);
	}
}

Graph GraphBuilder::toGraph(bool autoCompleteEdges, bool parallel) {
	Graph G(n, weighted, directed);
	if (name != "") {
		G.setName(name);
	}

	assert(G.outEdges.size() == n);
	assert(G.outEdgeWeights.size() == (weighted ? n : 0));
	assert(G.inEdges.size() == (directed ? n : 0));
	assert(G.inEdgeWeights.size() == ((weighted && directed) ? n : 0));

	// copy edges and weights
	if (autoCompleteEdges) {
		if (parallel) {
			toGraphParallel(G);
		} else {
			toGraphSequential(G);
		}
	} else {
		toGraphDirectSwap(G);
	}

	assert(G.outEdges.size() == n);
	assert(G.outEdgeWeights.size() == (weighted ? n : 0));
	assert(G.inEdges.size() == (directed ? n : 0));
	assert(G.inEdgeWeights.size() == ((weighted && directed) ? n : 0));

	G.m = numberOfEdges(G);
	G.shrinkToFit();

	reset();

	return G;
}

void GraphBuilder::toGraphDirectSwap(Graph &G) {
	G.outEdges = std::move(outEdges);
	G.outEdgeWeights = std::move(outEdgeWeights);
	G.inEdges = std::move(inEdges);
	G.inEdgeWeights = std::move(inEdgeWeights);
	if (!directed) {
		G.storedNumberOfSelfLoops = selfloops;
	} else if (selfloops % 2 == 0) {
		G.storedNumberOfSelfLoops = selfloops / 2;
	} else {
		throw std::runtime_error("Error, odd number of self loops added but each "
		                         "self loop must be added twice!");
	}
}

void GraphBuilder::toGraphParallel(Graph &G) {
	// basic idea of the parallelization:
	// 1) each threads collects its own data
	// 2) each node collects all its data from all threads

	int maxThreads = omp_get_max_threads();

	using Adjacencylists = std::vector<std::vector<storednode>>;
	using Weightlists = std::vector<std::vector<edgeweight>>;

	std::vector<Adjacencylists> inEdgesPerThread(maxThreads, Adjacencylists(n));
	std::vector<Weightlists> inWeightsPerThread(weighted ? maxThreads : 0,
	                                            Weightlists(n));
	std::vector<count> numberOfSelfLoopsPerThread(maxThreads, 0);

	// step 1
	parallelForNodes([&](node v) {
		int tid = omp_get_thread_num();
		for (index i = 0; i < outEdges[v].size(); i++) {
			node u = outEdges[v][i];
			if (directed || u != v) { // self loops don't need to be added twice in
				                        // undirected graphs
				inEdgesPerThread[tid][u].push_back(v);
				if (weighted) {
					edgeweight ew = outEdgeWeights[v][i];
					inWeightsPerThread[tid][u].push_back(ew);
				}
			}
			if (u == v) {
				numberOfSelfLoopsPerThread[tid]++;
			}
		}
	});

	// we have already half of the edges
	G.outEdges.swap(outEdges);
	G.outEdgeWeights.swap(outEdgeWeights);

	// step 2
	parallelForNodes([&](node v) {
		count inDeg = 0;
		count outDeg = G.outEdges[v].size();
		for (int tid = 0; tid < maxThreads; tid++) {
			inDeg += inEdgesPerThread[tid][v].size();
		}

		assert(inDeg <= n);
		assert(outDeg <= n);
		// allocate enough memory for all edges/weights
		if (directed) {
			G.inEdges[v].reserve(inDeg);
			if (weighted) {
				G.inEdgeWeights[v].reserve(inDeg);
			}
		} else {
			G.outEdges[v].reserve(outDeg + inDeg);
			if (weighted) {
				G.outEdgeWeights[v].reserve(outDeg + inDeg);
			}
		}

		// collect 'second' half of the edges
		if (directed) {
			for (int tid = 0; tid < maxThreads; tid++) {
				copyAndClear(inEdgesPerThread[tid][v], G.inEdges[v]);
			}
			if (weighted) {
				for (int tid = 0; tid < maxThreads; tid++) {
					copyAndClear(inWeightsPerThread[tid][v], G.inEdgeWeights[v]);
				}
			}

		} else {
			for (int tid = 0; tid < maxThreads; tid++) {
				copyAndClear(inEdgesPerThread[tid][v], G.outEdges[v]);
			}
			if (weighted) {
				for (int tid = 0; tid < maxThreads; tid++) {
					copyAndClear(inWeightsPerThread[tid][v], G.outEdgeWeights[v]);
				}
			}
		}
	});
	count numSelfLoops = 0;
#pragma omp parallel for reduction(+ : numSelfLoops)
	for (omp_index i = 0; i < static_cast<omp_index>(maxThreads); ++i)
		numSelfLoops += numberOfSelfLoopsPerThread[i];
	G.storedNumberOfSelfLoops = numSelfLoops;
}

void GraphBuilder::toGraphSequential(Graph &G) {
	std::vector<count> missingEdgesCounts(n, 0);
	count numberOfSelfLoops = 0;

	// 'first' half of the edges
	G.outEdges.swap(outEdges);
	G.outEdgeWeights.swap(outEdgeWeights);
	std::vector<count> outDeg(G.z);
	parallelForNodes([&](node v) { outDeg[v] = G.outEdges[v].size(); });

	// count missing edges for each node
	G.forNodes([&](node v) {
		// increase count of incoming edges for all neighbors
		for (node u : G.outEdges[v]) {
			if (directed || u != v) {
				++missingEdgesCounts[u];
			}
			if (u == v) {
				// self loops don't need to be added again
				// but we need to count them
				++numberOfSelfLoops;
			}
		}
	});

	// 'second' half the edges
	if (directed) {
		// directed: outEdges is complete, missing half edges are the inEdges
		// missingEdgesCounts are our inDegrees
		std::vector<count> inDeg = missingEdgesCounts;

		// reserve the exact amount of space needed
		G.forNodes([&](node v) {
			G.inEdges[v].reserve(inDeg[v]);
			if (weighted) {
				G.inEdgeWeights[v].reserve(inDeg[v]);
			}
		});

		// copy values
		G.forNodes([&](node v) {
			for (index i = 0; i < outDeg[v]; i++) {
				node u = G.outEdges[v][i];
				G.inEdges[u].push_back(v);
				if (weighted) {
					edgeweight ew = G.outEdgeWeights[v][i];
					G.inEdgeWeights[u].push_back(ew);
				}
			}
		});
	} else {
		// undirected: so far each edge is just saved at one node
		// add it to the other node as well

		// reserve the exact amount of space needed
		G.forNodes([&](node v) {
			G.outEdges[v].reserve(outDeg[v] + missingEdgesCounts[v]);
			if (weighted) {
				G.outEdgeWeights[v].reserve(outDeg[v] + missingEdgesCounts[v]);
			}
		});

		// copy values
		G.forNodes([&](node v) {
			// the first outDeg[v] edges in G.outEdges[v] are the first half edges
			// we are adding after outDeg[v]
			for (index i = 0; i < outDeg[v]; i++) {
				node u = G.outEdges[v][i];
				if (u != v) {
					G.outEdges[u].push_back(v);
					if (weighted) {
						edgeweight ew = G.outEdgeWeights[v][i];
						G.outEdgeWeights[u].push_back(ew);
					}
				} else {
					// ignore self loops here
				}
			}
		});
	}

	G.storedNumberOfSelfLoops = numberOfSelfLoops;
}

count GraphBuilder::numberOfEdges(const Graph &G) {
	count m = 0;
#pragma omp parallel for reduction(+ : m) if (n > (1 << 20))
	for (omp_index v = 0; v < static_cast<omp_index>(G.z); v++) {
		m += G.degree(v);
	}
	if (G.isDirected()) {
		return m;
	} else {
		// self loops are just counted once
		return (m - selfloops) / 2 + selfloops;
	}
}

} /* namespace NetworKit */
