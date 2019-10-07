/*
 * Graph.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter
 * (marvin.ritter@gmail.com)
 */

#include <cmath>
#include <map>
#include <random>
#include <sstream>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/** CONSTRUCTORS **/

Graph::Graph(count n, bool weighted, bool directed)
    : n(n), m(0), storedNumberOfSelfLoops(0), z(n), omega(0), t(0),

      weighted(weighted),  // indicates whether the graph is weighted or not
      directed(directed),  // indicates whether the graph is directed or not
      edgesIndexed(false), // edges are not indexed by default

      exists(n, true),

      /* for directed graphs inEdges stores an adjacencylist only considering
         incoming edges, for undirected graphs inEdges is not used*/
      inEdges(directed ? n : 0),

      /* for directed graphs outEdges stores an adjacencylist only considering
      outgoing edges, for undirected graphs outEdges stores the adjacencylist of
      undirected edges*/
      outEdges(n), inEdgeWeights(weighted && directed ? n : 0),
      outEdgeWeights(weighted ? n : 0), inEdgeIds(), outEdgeIds() {

	// set name from global id
	id = getNextGraphId();
	std::stringstream sstm;
	sstm << "G#" << id;
	name = sstm.str();
}

Graph::Graph(std::initializer_list<WeightedEdge> edges) : Graph(0, true) {
	using namespace std;

	/* Number of nodes = highest node index + 1 */
	for (const auto &edge : edges) {
		node x = max(edge.u, edge.v);
		while (numberOfNodes() <= x) {
			addNode();
		}
	}

	/* Now add all of the edges */
	for (const auto &edge : edges) {
		addEdge(edge.u, edge.v, edge.weight);
	}
}

Graph::Graph(const Graph &G, bool weighted, bool directed)
    : n(G.n), m(G.m), storedNumberOfSelfLoops(G.storedNumberOfSelfLoops),
      z(G.z), omega(0), t(G.t), weighted(weighted), directed(directed),
      edgesIndexed(false), // edges are not indexed by default
      exists(G.exists),

      // let the following be empty for the start, we fill them later
      inEdges(0), outEdges(0), inEdgeWeights(0),
      outEdgeWeights(0) {

	// set name from global id
	id = getNextGraphId();
	std::stringstream sstm;
	sstm << "G#" << id;
	name = sstm.str();

	if (G.isDirected() == directed) {
		// G.inEdges might be empty (if G is undirected), but
		// that's fine
		inEdges = G.inEdges;
		outEdges = G.outEdges;

		// copy weights if needed
		if (weighted) {
			if (G.isWeighted()) {
				// just copy from G, again either both graphs are directed or both are
				// undirected
				inEdgeWeights = G.inEdgeWeights;
				outEdgeWeights = G.outEdgeWeights;
			} else {
				// G has no weights, set defaultEdgeWeight for all edges
				if (directed) {
					inEdgeWeights.resize(z);
					for (node u = 0; u < z; u++) {
						inEdgeWeights[u] =
						    std::vector<edgeweight>(G.inEdges[u].size(), defaultEdgeWeight);
					}
				}

				outEdgeWeights.resize(z);
				for (node u = 0; u < z; u++) {
					outEdgeWeights[u] =
					    std::vector<edgeweight>(outEdges[u].size(), defaultEdgeWeight);
				}
			}
		}
	} else if (G.isDirected()) {
		// G is directed, but we want an undirected graph
		// so we need to combine the out and in stuff for every node
		outEdges.resize(z);
		for (node u = 0; u < z; u++) {

			// copy both out and in edges into our new outEdges
			outEdges[u].reserve(G.outEdges[u].size() + G.inEdges[u].size());
			outEdges[u].insert(outEdges[u].end(), G.outEdges[u].begin(),
			                   G.outEdges[u].end());
			outEdges[u].insert(outEdges[u].end(), G.inEdges[u].begin(),
			                   G.inEdges[u].end());
		}
		if (weighted) {
			if (G.isWeighted()) {
				// same for weights
				outEdgeWeights.resize(z);
				for (node u = 0; u < z; u++) {
					outEdgeWeights[u].reserve(G.outEdgeWeights[u].size() +
					                          G.inEdgeWeights[u].size());
					outEdgeWeights[u].insert(outEdgeWeights[u].end(),
					                         G.outEdgeWeights[u].begin(),
					                         G.outEdgeWeights[u].end());
					outEdgeWeights[u].insert(outEdgeWeights[u].end(),
					                         G.inEdgeWeights[u].begin(),
					                         G.inEdgeWeights[u].end());
				}
			} else {
				// we are undirected, so no need to write anything into inEdgeWeights
				outEdgeWeights.resize(z);
				for (node u = 0; u < z; u++) {
					outEdgeWeights[u] =
					    std::vector<edgeweight>(outEdges[u].size(), defaultEdgeWeight);
				}
			}
		}
	} else {
		// G is not directed, but this copy should be
		// generally we can can copy G.out stuff into our in stuff
		inEdges = G.outEdges;
		outEdges = G.outEdges;
		if (weighted) {
			if (G.isWeighted()) {
				inEdgeWeights = G.outEdgeWeights;
				outEdgeWeights = G.outEdgeWeights;
			} else {
				// initialize both inEdgeWeights and outEdgeWeights with the
				// defaultEdgeWeight
				inEdgeWeights.resize(z);
				for (node u = 0; u < z; u++) {
					inEdgeWeights[u] =
					    std::vector<edgeweight>(inEdges[u].size(), defaultEdgeWeight);
				}
				outEdgeWeights.resize(z);
				for (node u = 0; u < z; u++) {
					outEdgeWeights[u] =
					    std::vector<edgeweight>(outEdges[u].size(), defaultEdgeWeight);
				}
			}
		}
	}
}

void Graph::preallocateUndirected(node u, size_t size) {
	assert(!directed);
	assert(exists[u]);
	outEdges[u].reserve(size);
	if(weighted) {
		outEdgeWeights[u].reserve(size);
	}
	if(edgesIndexed) {
		outEdgeIds[u].reserve(size);
	}
}

void Graph::preallocateDirected(node u, size_t outSize, size_t inSize) {
	assert(directed);
	assert(exists[u]);
	inEdges[u].reserve(inSize);
	outEdges[u].reserve(outSize);

	if(weighted) {
		inEdgeWeights[u].reserve(inSize);
		outEdgeWeights[u].reserve(outSize);
	}
	if(edgesIndexed) {
		inEdgeIds[u].reserve(inSize);
		outEdgeIds[u].reserve(outSize);
	}
}
/** PRIVATE HELPERS **/

count Graph::getNextGraphId() {
	static count nextGraphId = 1;
	return nextGraphId++;
}

index Graph::indexInInEdgeArray(node v, node u) const {
	if (!directed) {
		return indexInOutEdgeArray(v, u);
	}
	for (index i = 0; i < inEdges[v].size(); i++) {
		node x = inEdges[v][i];
		if (x == u) {
			return i;
		}
	}
	return none;
}

index Graph::indexInOutEdgeArray(node u, node v) const {
	for (index i = 0; i < outEdges[u].size(); i++) {
		node x = outEdges[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

/** EDGE IDS **/

void Graph::indexEdges(bool force) {
	if (edgesIndexed && !force)
		return;

	omega = 0; // reset edge ids (for re-indexing)

	outEdgeIds.resize(outEdges.size());
	forNodes([&](node u) { outEdgeIds[u].resize(outEdges[u].size(), none); });

	if (directed) {
		inEdgeIds.resize(inEdges.size());
		forNodes([&](node u) { inEdgeIds[u].resize(inEdges[u].size(), none); });
	}

	// assign edge ids for edges in one direction
	forNodes([&](node u) {
		for (index i = 0; i < outEdges[u].size(); ++i) {
			node v = outEdges[u][i];
			if (v != none && (directed || (u >= v))) {
				// new id
				edgeid id = omega++;
				outEdgeIds[u][i] = id;
			}
		}
	});

	// copy edge ids for the edges in the other direction. Note that
	// "indexInOutEdgeArray" is slow which is why this second loop in parallel
	// makes sense.
	if (!directed) {
		balancedParallelForNodes([&](node u) {
			for (index i = 0; i < outEdges[u].size(); ++i) {
				node v = outEdges[u][i];
				if (v != none && outEdgeIds[u][i] == none) {
					index j = indexInOutEdgeArray(v, u);
					outEdgeIds[u][i] = outEdgeIds[v][j];
				}
			}
		});
	} else {
		balancedParallelForNodes([&](node u) {
			for (index i = 0; i < inEdges[u].size(); ++i) {
				node v = inEdges[u][i];
				if (v != none) {
					index j = indexInOutEdgeArray(v, u);
					inEdgeIds[u][i] = outEdgeIds[v][j];
				}
			}
		});
	}

	edgesIndexed = true; // remember that edges have been indexed so that addEdge
	                     // needs to create edge ids
}

edgeid Graph::edgeId(node u, node v) const {
	if (!edgesIndexed) {
		throw std::runtime_error(
		    "edges have not been indexed - call indexEdges first");
	}
	index i = indexInOutEdgeArray(u, v);

	if (i == none) {
		throw std::runtime_error("Edge does not exist");
	}

	edgeid id = outEdgeIds[u][i];
	return id;
}

/** GRAPH INFORMATION **/

std::string Graph::typ() const {
	if (weighted) {
		return directed ? "WeightedDirectedGraph" : "WeightedGraph";
	} else {
		return directed ? "DirectedGraph" : "Graph";
	}
}

void Graph::shrinkToFit() {
	exists.shrink_to_fit();

	inEdgeWeights.shrink_to_fit();
	for (auto &w : inEdgeWeights) {
		w.shrink_to_fit();
	}

	outEdgeWeights.shrink_to_fit();
	for (auto &w : outEdgeWeights) {
		w.shrink_to_fit();
	}

	inEdges.shrink_to_fit();
	for (auto &a : inEdges) {
		a.shrink_to_fit();
	}

	outEdges.shrink_to_fit();
	for (auto &a : outEdges) {
		a.shrink_to_fit();
	}
}

void Graph::compactEdges() {
	this->parallelForNodes([&](node u) {
		if (degreeOut(u) != outEdges[u].size()) {
			if (degreeOut(u) == 0) {
				outEdges[u].clear();
				if (weighted)
					outEdgeWeights[u].clear();
				if (edgesIndexed)
					outEdgeIds[u].clear();
			} else {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					while (i < outEdges[u].size() && outEdges[u][i] == none) {
						outEdges[u][i] = outEdges[u].back();
						outEdges[u].pop_back();

						if (weighted) {
							outEdgeWeights[u][i] = outEdgeWeights[u].back();
							outEdgeWeights[u].pop_back();
						}

						if (edgesIndexed) {
							outEdgeIds[u][i] = outEdgeIds[u].back();
							outEdgeIds[u].pop_back();
						}
					}
				}
			}
		}

		if (directed && degreeIn(u) != inEdges[u].size()) {
			if (degreeIn(u) == 0) {
				inEdges[u].clear();
				if (weighted)
					inEdgeWeights[u].clear();
				if (edgesIndexed)
					inEdgeIds[u].clear();
			} else {
				for (index i = 0; i < inEdges[u].size(); ++i) {
					while (i < inEdges[u].size() && inEdges[u][i] == none) {
						inEdges[u][i] = inEdges[u].back();
						inEdges[u].pop_back();

						if (weighted) {
							inEdgeWeights[u][i] = inEdgeWeights[u].back();
							inEdgeWeights[u].pop_back();
						}

						if (edgesIndexed) {
							inEdgeIds[u][i] = inEdgeIds[u].back();
							inEdgeIds[u].pop_back();
						}
					}
				}
			}
		}
	});
}

void Graph::sortEdges() {
	std::vector<std::vector<storednode>> targetAdjacencies(upperNodeIdBound());
	std::vector<std::vector<edgeweight>> targetWeight;
	std::vector<std::vector<edgeid>> targetEdgeIds;

	if (isWeighted()) {
		targetWeight.resize(upperNodeIdBound());
		forNodes([&](node u) { targetWeight[u].reserve(degree(u)); });
	}
	if (hasEdgeIds()) {
		targetEdgeIds.resize(upperNodeIdBound());
		forNodes([&](node u) { targetEdgeIds[u].reserve(degree(u)); });
	}

	forNodes([&](node u) { targetAdjacencies[u].reserve(degree(u)); });

	auto assignToTarget = [&](node u, node v, edgeweight w, edgeid eid) {
		targetAdjacencies[v].push_back(u);
		if (isWeighted()) {
			targetWeight[v].push_back(w);
		}
		if (hasEdgeIds()) {
			targetEdgeIds[v].push_back(eid);
		}
	};

	forNodes([&](node u) { forInEdgesOf(u, assignToTarget); });

	outEdges.swap(targetAdjacencies);
	outEdgeWeights.swap(targetWeight);
	outEdgeIds.swap(targetEdgeIds);

	if (isDirected()) {
		inEdges.swap(targetAdjacencies);
		inEdgeWeights.swap(targetWeight);
		inEdgeIds.swap(targetEdgeIds);

		forNodes([&](node u) {
			targetAdjacencies[u].resize(degreeIn(u));
			targetAdjacencies[u].shrink_to_fit();
			targetAdjacencies[u].clear();
			if (isWeighted()) {
				targetWeight[u].resize(degreeIn(u));
				targetWeight[u].shrink_to_fit();
				targetWeight[u].clear();
			}
			if (hasEdgeIds()) {
				targetEdgeIds[u].resize(degreeIn(u));
				targetEdgeIds[u].shrink_to_fit();
				targetEdgeIds[u].clear();
			}
		});

		forNodes([&](node u) { forEdgesOf(u, assignToTarget); });

		inEdges.swap(targetAdjacencies);
		inEdgeWeights.swap(targetWeight);
		inEdgeIds.swap(targetEdgeIds);
	}
}

count Graph::maxDegree() const { return computeMaxDegree(); }

count Graph::maxDegreeIn() const { return computeMaxDegree(true); }

count Graph::computeMaxDegree(const bool inDegree) const {
	count result = 0;
#ifndef NETWORKIT_OMP2
#pragma omp parallel for reduction(max : result)
	for (omp_index u = 0; u < upperNodeIdBound(); ++u) {
		result = std::max(result, inDegree ? degreeIn(u) : degreeOut(u));
	}
#else
	this->forNodes([&](const node u) {
		result = std::max(result, inDegree ? degreeIn(u) : degreeOut(u));
	});
#endif
	return result;
}

edgeweight Graph::maxWeightedDegree() const {
	return computeMaxWeightedDegree();
}

edgeweight Graph::maxWeightedDegreeIn() const {
	return computeMaxWeightedDegree(true);
}

edgeweight Graph::computeMaxWeightedDegree(const bool inDegree) const {
	edgeweight result = 0;
#ifndef NETWORKIT_OMP2
#pragma omp parallel for reduction(max : result)
	for (omp_index u = 0; u < upperNodeIdBound(); ++u) {
		result =
		    std::max(result, inDegree ? weightedDegreeIn(u) : weightedDegree(u));
	}
#else
	this->forNodes([&](const node u) {
		result =
		    std::max(result, inDegree ? weightedDegreeIn(u) : weightedDegree(u));
	});
#endif
	return result;
}

edgeweight Graph::computeWeightedDegree(const node &v,
                                        const bool inDegree) const {
	if (weighted) {
		edgeweight sum = 0.0;
		auto sumWeights = [&](node, edgeweight w) { sum += w; };
		if (inDegree) {
			forInNeighborsOf(v, sumWeights);
		} else {
			forNeighborsOf(v, sumWeights);
		}
		return sum;
	}
	return defaultEdgeWeight * (inDegree ? degreeIn(v) : degreeOut(v));
}

std::string Graph::toString() const {
	std::stringstream strm;
	strm << typ() << "(name=" << getName() << ", n=" << numberOfNodes()
	     << ", m=" << numberOfEdges() << ")";
	return strm.str();
}

/** COPYING **/

Graph Graph::copyNodes() const {
	Graph C(z, weighted, directed);
	for (node u = 0; u < z; ++u) {
		if (!exists[u]) {
			C.removeNode(u);
		}
	}
	return C;
}

/** NODE MODIFIERS **/

node Graph::addNode() {
	node v = z; // node gets maximum id
	z++;        // increment node range
	n++;        // increment node count
	assert(v == static_cast<node>(std::numeric_limits<storednode>::max()));

	// update per node data structures
	exists.push_back(true);

	// update per node data structures
	if (weighted) {

		std::vector<edgeweight> edgeWeight;
		inEdgeWeights.push_back(edgeWeight);
		outEdgeWeights.push_back(edgeWeight);
	}

	outEdges.push_back(std::vector<storednode>{});
	if (directed) {
		inEdges.push_back(std::vector<storednode>{});
	}

	return v;
}

node Graph::addNode(float x, float y) {
	node v = addNode();
	std::vector<float> coords = {x, y};
	coordinates.addCoordinates(coords);
	return v;
}

void Graph::removeNode(node v) {
	assert(v < z);
	assert(exists[v]);

	// Remove all outgoing and ingoing edges
	while (!outEdges[v].empty())
		removeEdge(v, outEdges[v].front());
	if (isDirected())
		while (!inEdges[v].empty())
			removeEdge(inEdges[v].front(), v);

	exists[v] = false;
	n--;
}

void Graph::restoreNode(node v) {
	assert(v < z);
	assert(!exists[v]);

	exists[v] = true;
	n++;
}

/** NODE PROPERTIES **/

edgeweight Graph::weightedDegree(const node &v) const {
	return computeWeightedDegree(v);
}

edgeweight Graph::weightedDegreeIn(const node &v) const {
	return computeWeightedDegree(v, true);
}

edgeweight Graph::volume(node v) const {
	if (weighted) {
		edgeweight sum = 0.0;
		for (index i = 0; i < outEdges[v].size(); i++) {
			node u = outEdges[v][i];
			if (u == v) {
				sum += 2 * outEdgeWeights[v][i];
			} else if (u != none) {
				sum += outEdgeWeights[v][i];
			}
		}
		return sum;
	} else {
		count c = outEdges[v].size();
		for (node u : outEdges[v]) {
			if (u == v) {
				c++;
			}
		}
		return c * defaultEdgeWeight;
	}
}

node Graph::randomNode() const {
	if (numberOfNodes() == 0) {
		return none;
	}

	node v;
	do {
		v = Aux::Random::integer(z - 1);
	} while (!exists[v]);

	return v;
}

node Graph::randomNeighbor(node u) const {
	if (outEdges[u].empty()) {
		return none;
	}

	node v;
	do {
		index i = Aux::Random::integer(outEdges[u].size() - 1);
		v = outEdges[u][i];
	} while (v == none);

	return v;
}

/** EDGE MODIFIERS **/

void Graph::addEdge(node u, node v, edgeweight ew) {
	assert(u < z);
	assert(exists[u]);
	assert(v < z);
	assert(exists[v]);

	// increase number of edges
	++m;
	outEdges[u].push_back(v);

	// if edges indexed, give new id
	if (edgesIndexed) {
		edgeid id = omega++;
		outEdgeIds[u].push_back(id);
	}

	if (directed) {
		inEdges[v].push_back(u);

		if (edgesIndexed) {
			inEdgeIds[v].push_back(id);
		}

		if (weighted) {
			inEdgeWeights[v].push_back(ew);
			outEdgeWeights[u].push_back(ew);
		}

	} else if (u == v) { // self-loop case
		if (weighted) {
			outEdgeWeights[u].push_back(ew);
		}
	} else { // undirected, no self-loop
		outEdges[v].push_back(u);

		if (weighted) {
			outEdgeWeights[u].push_back(ew);
			outEdgeWeights[v].push_back(ew);
		}

		if (edgesIndexed) {
			outEdgeIds[v].push_back(omega - 1);
		}
	}

	if (u == v) { // count self loop
		++storedNumberOfSelfLoops;
	}
}
void Graph::addPartialEdge(Unsafe, node u, node v, edgeweight ew, uint64_t index) {
	assert(u < z);
	assert(exists[u]);
	assert(v < z);
	assert(exists[v]);

	outEdges[u].push_back(v);

	// if edges indexed, give new id
	if (edgesIndexed) {
		outEdgeIds[u].push_back(index);
	}
	if (weighted) {
		outEdgeWeights[u].push_back(ew);
	}
}
void Graph::addPartialOutEdge(Unsafe, node u, node v, edgeweight ew, uint64_t index) {
	assert(u < z);
	assert(exists[u]);
	assert(v < z);
	assert(exists[v]);

	outEdges[u].push_back(v);

	// if edges indexed, give new id
	if (edgesIndexed) {
		outEdgeIds[u].push_back(index);
	}
	if (weighted) {
		outEdgeWeights[u].push_back(ew);
	}
}
void Graph::addPartialInEdge(Unsafe, node u, node v, edgeweight ew, uint64_t index) {
	assert(u < z);
	assert(exists[u]);
	assert(v < z);
	assert(exists[v]);

	inEdges[u].push_back(v);

	if (edgesIndexed) {
		inEdgeIds[u].push_back(index);
	}
	if (weighted) {
		inEdgeWeights[u].push_back(ew);
	}
}

template <typename T>
void erase(node u, index idx, std::vector<std::vector<T>> &vec) {
	vec[u][idx] = vec[u].back();
	vec[u].pop_back();
}

void Graph::removeEdge(node u, node v) {
	assert(u < z);
	assert(exists[u]);
	assert(v < z);
	assert(exists[v]);

	index vi = indexInOutEdgeArray(u, v);
	index ui = indexInInEdgeArray(v, u);

	if (vi == none) {
		std::stringstream strm;
		strm << "edge (" << u << "," << v << ") does not exist";
		throw std::runtime_error(strm.str());
	}

	m--; // decrease number of edges
	erase<storednode>(u, vi, outEdges);
	if (weighted) {
		erase<edgeweight>(u, vi, outEdgeWeights);
	}

	if (directed) {
		assert(ui != none);

		erase<storednode>(v, ui, inEdges);
		if (weighted) {
			erase<edgeweight>(v, ui, inEdgeWeights);
		}
	} else if (u != v) {
		// undirected, not self-loop
		erase<storednode>(v, ui, outEdges);
		if (weighted) {
			erase<edgeweight>(v, ui, outEdgeWeights);
		}
	}

	if (u == v) {
		storedNumberOfSelfLoops--;
	}

	// dose not make a lot of sense do remove attributes,
	// cause the edge is marked as deleted and we have no null values for the
	// attributes
}

void Graph::removeAllEdges() {
#pragma omp parallel for
	for (omp_index u = 0; u < z; ++u) {
		outEdges[u].clear();
		if (isWeighted()) {
			outEdgeWeights[u].clear();
		}
		if (isDirected()) {
			inEdges[u].clear();
			if (isWeighted()) {
				inEdgeWeights[u].clear();
			}
		}
	}

	m = 0;
}

void Graph::removeEdgesFromIsolatedSet(const std::vector<node> &nodesInSet) {
	count removedEdges = 0;
	for (node u : nodesInSet) {
		removedEdges += outEdges[u].size();
		outEdges[u].clear();
		if (this->weighted) {
			outEdgeWeights[u].clear();
		}
		if (this->isDirected()) {
			inEdges[u].clear();
		}
	}
	this->m -= removedEdges;
}

void Graph::removeSelfLoops() {
#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(z); ++i) {
		const auto u = static_cast<node>(i);
		auto isSelfLoop = [u](node v) { return u == v; };
		removeAdjacentEdges(u, isSelfLoop);
		if (isDirected()) {
			removeAdjacentEdges(u, isSelfLoop, true);
		}
	}

	m -= storedNumberOfSelfLoops;
	storedNumberOfSelfLoops = 0;
}

void Graph::removeMultiEdges() {
	count removedEdges = 0;
	count removedSelfLoops = 0;
	std::unordered_set<node> nodes;

	forNodes([&](node u) {
		nodes.reserve(degree(u));
		auto isMultiedge = [&nodes](node v) { return !nodes.insert(v).second; };
		auto result = removeAdjacentEdges(u, isMultiedge);
		removedEdges += result.first;
		removedSelfLoops += result.second;
		if (isDirected()) {
			nodes.clear();
			removeAdjacentEdges(u, isMultiedge, true).first;
		}
		nodes.clear();
	});

	if (!isDirected()) {
		assert(!(removedEdges % 2));
		removedEdges /= 2;
	}

	m -= removedEdges + removedSelfLoops;
	storedNumberOfSelfLoops -= removedSelfLoops;
}

void Graph::swapEdge(node s1, node t1, node s2, node t2) {
	index s1t1 = indexInOutEdgeArray(s1, t1);
	if (s1t1 == none)
		throw std::runtime_error("The first edge does not exist");
	index t1s1 = indexInInEdgeArray(t1, s1);

	index s2t2 = indexInOutEdgeArray(s2, t2);
	if (s2t2 == none)
		throw std::runtime_error("The second edge does not exist");
	index t2s2 = indexInInEdgeArray(t2, s2);

	std::swap(outEdges[s1][s1t1], outEdges[s2][s2t2]);

	if (directed) {
		std::swap(inEdges[t1][t1s1], inEdges[t2][t2s2]);

		if (weighted) {
			std::swap(inEdgeWeights[t1][t1s1], inEdgeWeights[t2][t2s2]);
		}

		if (edgesIndexed) {
			std::swap(inEdgeIds[t1][t1s1], inEdgeIds[t2][t2s2]);
		}
	} else {
		std::swap(outEdges[t1][t1s1], outEdges[t2][t2s2]);

		if (weighted) {
			std::swap(outEdgeWeights[t1][t1s1], outEdgeWeights[t2][t2s2]);
		}

		if (edgesIndexed) {
			std::swap(outEdgeIds[t1][t1s1], outEdgeIds[t2][t2s2]);
		}
	}
}

bool Graph::hasEdge(node u, node v) const {
	if (!directed &&
	    outEdges[u].size() > outEdges[v].size()) {
		return indexInOutEdgeArray(v, u) != none;
	} else if (directed &&
	           outEdges[u].size() > inEdges[v].size()) {
		return indexInInEdgeArray(v, u) != none;
	} else {
		return indexInOutEdgeArray(u, v) != none;
	}
}

std::pair<node, node> Graph::randomEdge(bool uniformDistribution) const {
	if (m == 0) {
		throw std::runtime_error("Error: the graph has ne edges!");
	}

	if (uniformDistribution) {
		/*
		 * The simple idea here is to interpret all neighborhoods next to each other, resulting
		 * in a virtual vector of size m. Then we draw a random index and return the edge.
		 * For undirected edges, the vector has size 2m; but the idea remains. There is one minor
		 * complication for undirected edges with self-loops: each edge {u,v} with u != v is stored
		 * twice (once in the neighborhood of u, once in v) but a loop (u, u) is only stored once.
		 * To equalize the probabilities we reject edges {u,v} with u > v and try again. This leads
		 * to less than two expected trails in and is only done for undirected graphs with self-loops.
		 */

		while (true) {
			const auto upper = directed
				? numberOfEdges()
				: 2 * numberOfEdges() - numberOfSelfLoops();
			auto idx = Aux::Random::index(upper);

			node u, v;

			if (idx > upper / 2) {
				// assuming degrees are somewhat distributed uniformly, it's better to start with
				// larger nodes for large indices. In this case we have to mirror the index:
				idx  = (upper-1) - idx;

				for(u = upperNodeIdBound() - 1; idx >= degreeOut(u); --u) {
					idx -= degreeOut(u);
				}

				v = outEdges[u][outEdges[u].size() - 1 - idx];

			} else {
				for(u = 0; idx >= degreeOut(u); ++u) {
					assert(u < upperNodeIdBound());
					idx -= degreeOut(u);
				}

				v = outEdges[u][idx];

			}

			if (numberOfSelfLoops() && !directed && u > v)
				// reject (see above)
				continue;

			return {u, v};
		}
	}

	node u; // we will return edge (u, v)

	// fast way, but not a uniform random edge!
	do {
		u = randomNode();
	} while (outEdges[u].empty());

	const auto v = randomNeighbor(u);

	return {u, v};
}

std::vector<std::pair<node, node>> Graph::randomEdges(count nr) const {
	if (numberOfEdges() == 0) {
		throw std::runtime_error(
		    "Graph has no edges to sample from. Add edges to the graph first.");
	}
	std::vector<std::pair<node, node>> edges;

	auto& gen = Aux::Random::getURNG();
	std::vector<count> outDeg(upperNodeIdBound());
	for (count i = 0; i < upperNodeIdBound(); ++i) {
		outDeg[i] = outEdges[i].size();
	}
	std::discrete_distribution<count> distribution(outDeg.begin(), outDeg.end());

	for (index i = 0; i < nr; i++) {
		node u, v; // we will pick edge (u, v)
		if (directed) {
			u = distribution(gen);
			// should always be the case as  without
			// edges should have probability 0
			assert(outEdges[u].size() > 0);
			v = randomNeighbor(u);
		} else {
			// self-loops which appear only once in the outEdge arrays
			// easiest way it to ignore edges (u, v) with u > v
			do {
				u = distribution(gen);
				// should always be the case as  without
				// edges should have probability 0
				assert(outEdges[u].size() > 0);
				v = randomNeighbor(u);
			} while (u > v);
		}
		edges.push_back({u, v});
	}

	return edges;
}

/** GLOBAL PROPERTIES **/

count Graph::numberOfSelfLoops() const { return storedNumberOfSelfLoops; }

/** EDGE ATTRIBUTES **/

edgeweight Graph::weight(node u, node v) const {
	index vi = indexInOutEdgeArray(u, v);
	if (vi == none) {
		return nullWeight;
	} else {
		return weighted ? outEdgeWeights[u][vi] : defaultEdgeWeight;
	}
}

void Graph::setWeight(node u, node v, edgeweight ew) {
	if (!weighted) {
		throw std::runtime_error("Cannot set edge weight in unweighted graph.");
	}

	index vi = indexInOutEdgeArray(u, v);
	if (vi == none) {
		// edge does not exist, create it, but warn user
		TRACE("Setting edge weight of a nonexisting edge will create the edge.");
		addEdge(u, v, ew);
		return;
	}

	outEdgeWeights[u][vi] = ew;
	if (directed) {
		index ui = indexInInEdgeArray(v, u);
		inEdgeWeights[v][ui] = ew;
	} else if (u != v) {
		index ui = indexInInEdgeArray(v, u);
		outEdgeWeights[v][ui] = ew;
	}
}

void Graph::increaseWeight(node u, node v, edgeweight ew) {
	if (!weighted) {
		throw std::runtime_error(
		    "Cannot increase edge weight in unweighted graph.");
	}

	index vi = indexInOutEdgeArray(u, v);
	if (vi == none) {
		// edge does not exits, create it, but warn user
		addEdge(u, v, ew);
		return;
	}

	outEdgeWeights[u][vi] += ew;
	if (directed) {
		index ui = indexInInEdgeArray(v, u);
		inEdgeWeights[v][ui] += ew;
	} else if (u != v) {
		index ui = indexInInEdgeArray(v, u);
		outEdgeWeights[v][ui] += ew;
	}
}

/** SUMS **/

edgeweight Graph::totalEdgeWeight() const {
	if (weighted) {
		edgeweight sum = 0.0;
		forEdges([&](node, node, edgeweight ew) { sum += ew; });
		return sum;
	} else {
		return numberOfEdges() * defaultEdgeWeight;
	}
}

/** Collections **/

std::vector<node> Graph::nodes() const {
	std::vector<node> nodes;
	nodes.reserve(numberOfNodes());
	this->forNodes([&](node u) { nodes.push_back(u); });
	return nodes;
}

std::vector<std::pair<node, node>> Graph::edges() const {
	std::vector<std::pair<node, node>> edges;
	edges.reserve(numberOfEdges());
	this->forEdges(
	    [&](node u, node v) { edges.push_back(std::pair<node, node>(u, v)); });
	return edges;
}

std::vector<node> Graph::neighbors(node u) const {
	std::vector<node> neighbors;
	neighbors.reserve(degree(u));
	this->forNeighborsOf(u, [&](node v) { neighbors.push_back(v); });
	return neighbors;
}

Graph Graph::transpose() const {
	if (directed == false) {
		throw std::runtime_error("The transpose of an undirected graph is "
								 "identical to the original graph.");
	}

	Graph GTranspose(z, weighted, true);

	// prepare edge id storage if input has indexed edges
	GTranspose.edgesIndexed = edgesIndexed;
	GTranspose.omega = omega;
	if (edgesIndexed) {
		GTranspose.inEdgeIds.resize(z);
		GTranspose.outEdgeIds.resize(z);
	}

	#pragma omp parallel for
	for (omp_index u = 0; u <  static_cast<omp_index>(z); ++u) {
		if (exists[u]) {
			GTranspose.preallocateDirected(u, degreeIn(u), degreeOut(u));

			forInEdgesOf(u, [&] (node, node v, edgeweight w, edgeid id) {
				GTranspose.addPartialOutEdge(unsafe, u, v, w, id);
			});

			forEdgesOf(u, [&] (node, node v, edgeweight w, edgeid id) {
				GTranspose.addPartialInEdge(unsafe, u, v, w, id);
			});

		} else {
			#pragma omp critical
			GTranspose.removeNode(u);
		}
	}

	GTranspose.setEdgeCount(unsafe, numberOfEdges());
	GTranspose.setNumberOfSelfLoops(unsafe, numberOfSelfLoops());
	GTranspose.t = t;
	GTranspose.setName(getName() + "Transpose");

	assert(GTranspose.checkConsistency());

	return GTranspose;
}

Graph Graph::toUndirected() const {
	if (directed == false) {
		throw std::runtime_error("this graph is already undirected");
	}
	Graph U(*this, weighted, false);
	return U;
}

Graph Graph::toUnweighted() const {
	if (weighted == false) {
		throw std::runtime_error("this graph is already unweighted");
	}
	Graph U(*this, false, directed);
	return U;
}

bool Graph::checkConsistency() const {
	// check for multi-edges
	std::vector<node> lastSeen(z, none);
	bool noMultiEdges = true;
	auto noMultiEdgesDetected = [&noMultiEdges]() { return noMultiEdges; };
	forNodesWhile(noMultiEdgesDetected, [&](node v) {
		forNeighborsOf(v, [&](node u) {
			if (lastSeen[u] == v) {
				noMultiEdges = false;
				DEBUG("Multiedge found between ", u, " and ", v, "!");
			}
			lastSeen[u] = v;
		});
	});

	return noMultiEdges;
}

void Graph::append(const Graph &G) {
	std::map<node, node> nodeMap;
	G.forNodes([&](node u) {
		node u_ = this->addNode();
		nodeMap[u] = u_;
	});
	if (this->isWeighted()) {
		G.forEdges([&](node u, node v, edgeweight ew) {
			this->addEdge(nodeMap[u], nodeMap[v], ew);
		});
	} else {
		G.forEdges([&](node u, node v) { this->addEdge(nodeMap[u], nodeMap[v]); });
	}
}

void Graph::merge(const Graph &G) {
	// TODO: handle edge weights
	G.forEdges([&](node u, node v) {
		// naive implementation takes $O(m \cdot d)$ for $m$ edges and max. degree
		// $d$ in this graph
		if (!this->hasEdge(u, v)) {
			this->addEdge(u, v);
		}
	});
}

// SUBGRAPHS

Graph Graph::subgraphFromNodes(const std::unordered_set<node> &nodes, bool includeOutNeighbors, bool includeInNeighbors) const {
	const auto neighbors = [&] {
		std::unordered_set<node> neighbors;

		if (!includeOutNeighbors && !includeInNeighbors)
			return neighbors;

		for (node u : nodes) {
			if (includeOutNeighbors)
				for(const node v : neighborRange(u))
					neighbors.insert(v);

			if (includeInNeighbors)
				for(const node v : inNeighborRange(u))
					neighbors.insert(v);
		}

		return neighbors;
	}();

	/*
	 * returns one of three different relevance scores:
	 * 2: Is in the original nodes set
	 * 1: Is a relevant neighbor (i.e., respective include*Neighbor was set)
	 * 0: Neither of both
	 */
	auto isRelevantNode = [&] (const node u) {
		if (nodes.find(u) != nodes.end()) return 2;
		if (!neighbors.empty() && neighbors.find(u) != neighbors.end()) return 1;
		return 0;
	};

	Graph S(upperNodeIdBound(), isWeighted(), isDirected());
	// delete all nodes that are not in the node set
	S.forNodes([&](node u) {
		if (!isRelevantNode(u)) {
			S.removeNode(u);
		}
	});

	forEdges([&](node u, node v, edgeweight w) {
		// only include edges if at least one endpoint is in nodes (relevance 2),
		// and the other is either in nodes or in neighbors (relevance >= 1)
		if (isRelevantNode(u) + isRelevantNode(v) > 2) {
			S.addEdge(u, v, w);
		}
	});

	return S;
}

} /* namespace NetworKit */
