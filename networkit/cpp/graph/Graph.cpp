/*
 * Graph.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Christian Staudt
 *              Klara Reichard <klara.reichard@gmail.com>
 *              Marvin Ritter <marvin.ritter@gmail.com>
 */

// networkit-format

#include <cmath>
#include <map>
#include <random>
#include <sstream>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

/** CONSTRUCTORS **/

Graph::Graph(count n, bool weighted, bool directed)
    : n(n), m(0), storedNumberOfSelfLoops(0), z(n), omega(0), t(0),

      weighted(weighted),  // indicates whether the graph is weighted or not
      directed(directed),  // indicates whether the graph is directed or not
      edgesIndexed(false), // edges are not indexed by default

      exists(n, true),

      /* for directed graphs inEdges stores an adjacency list only considering
         incoming edges, for undirected graphs inEdges is not used*/
      inEdges(directed ? n : 0),

      /* for directed graphs outEdges stores an adjacency list only considering
      outgoing edges, for undirected graphs outEdges stores the adjacency list of
      undirected edges*/
      outEdges(n), inEdgeWeights(weighted && directed ? n : 0), outEdgeWeights(weighted ? n : 0),
      inEdgeIds(), outEdgeIds() {}

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
    : n(G.n), m(G.m), storedNumberOfSelfLoops(G.storedNumberOfSelfLoops), z(G.z), omega(0), t(G.t),
      weighted(weighted), directed(directed),
      edgesIndexed(false), // edges are not indexed by default
      exists(G.exists),

      // let the following be empty for the start, we fill them later
      inEdges(0), outEdges(0), inEdgeWeights(0), outEdgeWeights(0) {

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
            outEdges[u].insert(outEdges[u].end(), G.outEdges[u].begin(), G.outEdges[u].end());
            outEdges[u].insert(outEdges[u].end(), G.inEdges[u].begin(), G.inEdges[u].end());
        }
        if (weighted) {
            if (G.isWeighted()) {
                // same for weights
                outEdgeWeights.resize(z);
                for (node u = 0; u < z; u++) {
                    outEdgeWeights[u].reserve(G.outEdgeWeights[u].size()
                                              + G.inEdgeWeights[u].size());
                    outEdgeWeights[u].insert(outEdgeWeights[u].end(), G.outEdgeWeights[u].begin(),
                                             G.outEdgeWeights[u].end());
                    outEdgeWeights[u].insert(outEdgeWeights[u].end(), G.inEdgeWeights[u].begin(),
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
    if (weighted) {
        outEdgeWeights[u].reserve(size);
    }
    if (edgesIndexed) {
        outEdgeIds[u].reserve(size);
    }
}

void Graph::preallocateDirected(node u, size_t outSize, size_t inSize) {
    preallocateDirectedOutEdges(u, outSize);
    preallocateDirectedInEdges(u, inSize);
}

void Graph::preallocateDirectedOutEdges(node u, size_t outSize) {
    assert(directed);
    assert(exists[u]);
    outEdges[u].reserve(outSize);

    if (weighted) {
        outEdgeWeights[u].reserve(outSize);
    }
    if (edgesIndexed) {
        outEdges[u].reserve(outSize);
    }
}

void Graph::preallocateDirectedInEdges(node u, size_t inSize) {
    assert(directed);
    assert(exists[u]);
    inEdges[u].reserve(inSize);

    if (weighted) {
        inEdgeWeights[u].reserve(inSize);
    }
    if (edgesIndexed) {
        inEdgeIds[u].reserve(inSize);
    }
}
/** PRIVATE HELPERS **/

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
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    index i = indexInOutEdgeArray(u, v);

    if (i == none) {
        throw std::runtime_error("Edge does not exist");
    }
    edgeid id = outEdgeIds[u][i];
    return id;
}

/** GRAPH INFORMATION **/

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
    std::vector<std::vector<node>> targetAdjacencies(upperNodeIdBound());
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

edgeweight Graph::computeWeightedDegree(node u, bool inDegree, bool countSelfLoopsTwice) const {
    if (weighted) {
        edgeweight sum = 0.0;
        auto sumWeights = [&](node v, edgeweight w) {
            sum += (countSelfLoopsTwice && u == v) ? 2. * w : w;
        };
        if (inDegree) {
            forInNeighborsOf(u, sumWeights);
        } else {
            forNeighborsOf(u, sumWeights);
        }
        return sum;
    }

    count sum = inDegree ? degreeIn(u) : degreeOut(u);
    auto countSelfLoops = [&](node v) { sum += (u == v); };

    if (countSelfLoopsTwice && numberOfSelfLoops()) {
        if (inDegree) {
            forInNeighborsOf(u, countSelfLoops);
        } else {
            forNeighborsOf(u, countSelfLoops);
        }
    }

    return static_cast<edgeweight>(sum);
}

/** NODE MODIFIERS **/

node Graph::addNode() {
    node v = z; // node gets maximum id
    z++;        // increment node range
    n++;        // increment node count

    // update per node data structures
    exists.push_back(true);

    outEdges.emplace_back();
    if (weighted)
        outEdgeWeights.emplace_back();
    if (edgesIndexed)
        outEdgeIds.emplace_back();

    if (directed) {
        inEdges.emplace_back();
        if (weighted)
            inEdgeWeights.emplace_back();
        if (edgesIndexed)
            inEdgeIds.emplace_back();
    }

    return v;
}

node Graph::addNodes(count numberOfNewNodes) {
    if (numberOfNewNodes < 10) {
        // benchmarks suggested, it's cheaper to call 10 time emplace_back than resizing.
        while (numberOfNewNodes--)
            addNode();

        return z - 1;
    }

    z += numberOfNewNodes;
    n += numberOfNewNodes;

    // update per node data structures
    exists.resize(z, true);

    outEdges.resize(z);
    if (weighted)
        outEdgeWeights.resize(z);
    if (edgesIndexed)
        outEdgeIds.resize(z);

    if (directed) {
        inEdges.resize(z);
        if (weighted)
            inEdgeWeights.resize(z);
        if (edgesIndexed)
            inEdgeIds.resize(z);
    }

    return z - 1;
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

edgeweight Graph::weightedDegree(node u, bool countSelfLoopsTwice) const {
    return computeWeightedDegree(u, false, countSelfLoopsTwice);
}

edgeweight Graph::weightedDegreeIn(node u, bool countSelfLoopsTwice) const {
    return computeWeightedDegree(u, true, countSelfLoopsTwice);
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
            inEdgeIds[v].push_back(omega - 1);
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

    const auto isLoop = (u == v);
    m--; // decrease number of edges
    if (isLoop)
        storedNumberOfSelfLoops--;

    erase<node>(u, vi, outEdges);
    if (weighted) {
        erase<edgeweight>(u, vi, outEdgeWeights);
    }
    if (edgesIndexed) {
        erase<edgeid>(u, vi, outEdgeIds);
    }

    if (directed) {
        assert(ui != none);

        erase<node>(v, ui, inEdges);
        if (weighted) {
            erase<edgeweight>(v, ui, inEdgeWeights);
        }
        if (edgesIndexed) {
            erase<edgeid>(v, ui, inEdgeIds);
        }
    } else if (!isLoop) {
        // undirected, not self-loop
        erase<node>(v, ui, outEdges);
        if (weighted) {
            erase<edgeweight>(v, ui, outEdgeWeights);
        }
        if (edgesIndexed) {
            erase<edgeid>(v, ui, outEdgeIds);
        }
    }
}

void Graph::removeAllEdges() {
    parallelForNodes([&](const node u) {
        removePartialOutEdges(unsafe, u);
        if (isDirected()) {
            removePartialInEdges(unsafe, u);
        }
    });

    m = 0;
}

void Graph::removeSelfLoops() {
    parallelForNodes([&](const node u) {
        auto isSelfLoop = [u](const node v) { return u == v; };
        removeAdjacentEdges(u, isSelfLoop);
        if (isDirected()) {
            removeAdjacentEdges(u, isSelfLoop, true);
        }
    });

    m -= storedNumberOfSelfLoops;
    storedNumberOfSelfLoops = 0;
}

void Graph::removeMultiEdges() {
    count removedEdges = 0;
    count removedSelfLoops = 0;
    std::unordered_set<node> nodes;

    forNodes([&](const node u) {
        nodes.reserve(degree(u));
        auto isMultiedge = [&nodes](const node v) { return !nodes.insert(v).second; };
        auto result = removeAdjacentEdges(u, isMultiedge);
        removedEdges += result.first;
        removedSelfLoops += result.second;
        if (isDirected()) {
            nodes.clear();
            removeAdjacentEdges(u, isMultiedge, true);
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

bool Graph::hasEdge(node u, node v) const noexcept {
    if (u >= z || v >= z) {
        return false;
    }
    if (!directed && outEdges[u].size() > outEdges[v].size()) {
        return indexInOutEdgeArray(v, u) != none;
    } else if (directed && outEdges[u].size() > inEdges[v].size()) {
        return indexInInEdgeArray(v, u) != none;
    } else {
        return indexInOutEdgeArray(u, v) != none;
    }
}

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
        throw std::runtime_error("Cannot increase edge weight in unweighted graph.");
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

void Graph::setWeightAtIthNeighbor(Unsafe, node u, index i, edgeweight ew) {
    outEdgeWeights[u][i] = ew;
}

/** SUMS **/

edgeweight Graph::totalEdgeWeight() const noexcept {
    if (weighted) {
        edgeweight sum = 0.0;
        forEdges([&](node, node, edgeweight ew) { sum += ew; });
        return sum;
    } else {
        return numberOfEdges() * defaultEdgeWeight;
    }
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

    bool correctNodeUpperbound = (z == outEdges.size()) && ((directed ? z : 0) == inEdges.size())
                                 && ((weighted ? z : 0) == outEdgeWeights.size())
                                 && ((weighted && directed ? z : 0) == inEdgeWeights.size())
                                 && ((edgesIndexed ? z : 0) == outEdgeIds.size())
                                 && ((edgesIndexed && directed ? z : 0) == inEdgeIds.size());

    if (!correctNodeUpperbound)
        DEBUG("Saved node upper bound doesn't actually match the actual node upper bound!");

    count NumberOfOutEdges = 0;
    count NumberOfOutEdgeWeights = 0;
    count NumberOfOutEdgeIds = 0;
    for (index i = 0; i < outEdges.size(); i++) {
        NumberOfOutEdges += outEdges[i].size();
    }
    if (weighted)
        for (index i = 0; i < outEdgeWeights.size(); i++) {
            NumberOfOutEdgeWeights += outEdgeWeights[i].size();
        }
    if (edgesIndexed)
        for (index i = 0; i < outEdgeIds.size(); i++) {
            NumberOfOutEdgeIds += outEdgeIds[i].size();
        }

    count NumberOfInEdges = 0;
    count NumberOfInEdgeWeights = 0;
    count NumberOfInEdgeIds = 0;
    if (directed) {
        for (index i = 0; i < inEdges.size(); i++) {
            NumberOfInEdges += inEdges[i].size();
        }
        if (weighted)
            for (index i = 0; i < inEdgeWeights.size(); i++) {
                NumberOfInEdgeWeights += inEdgeWeights[i].size();
            }
        if (edgesIndexed)
            for (index i = 0; i < inEdgeIds.size(); i++) {
                NumberOfInEdgeIds += inEdgeIds[i].size();
            }
    }

    if (!directed) {
        NumberOfOutEdges = (NumberOfOutEdges + storedNumberOfSelfLoops) / 2;
        if (weighted)
            NumberOfOutEdgeWeights = (NumberOfOutEdgeWeights + storedNumberOfSelfLoops) / 2;
        if (edgesIndexed)
            NumberOfOutEdgeIds = (NumberOfOutEdgeIds + storedNumberOfSelfLoops) / 2;
    }

    bool correctNumberOfEdges = (m == NumberOfOutEdges) && ((directed ? m : 0) == NumberOfInEdges)
                                && ((weighted ? m : 0) == NumberOfOutEdgeWeights)
                                && ((weighted && directed ? m : 0) == NumberOfInEdgeWeights)
                                && ((edgesIndexed ? m : 0) == NumberOfOutEdgeIds)
                                && ((edgesIndexed && directed ? m : 0) == NumberOfInEdgeIds);

    if (!correctNumberOfEdges)
        DEBUG("Saved number of edges is incorrect!");

    return noMultiEdges && correctNodeUpperbound && correctNumberOfEdges;
}

} /* namespace NetworKit */
