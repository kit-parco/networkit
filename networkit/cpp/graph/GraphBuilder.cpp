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

GraphBuilder::GraphBuilder(count n, bool weighted, bool directed)
    : n(n), selfloops(0), weighted(weighted), directed(directed),
      outEdges(n), outEdgeWeights(weighted ? n : 0), inEdges(directed ? n : 0),
      inEdgeWeights((directed && weighted) ? n : 0) {}

void GraphBuilder::reset(count n) {
    this->n = n;
    selfloops = 0;
    outEdges.assign(n, std::vector<node>{});
    outEdgeWeights.assign(isWeighted() ? n : 0, std::vector<edgeweight>{}),
        inEdges.assign(isDirected() ? n : 0, std::vector<node>{}),
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
    outEdges.emplace_back();
    if (weighted) {
        outEdgeWeights.emplace_back();
    }
    if (directed) {
        inEdges.emplace_back();
        if (weighted) {
            inEdgeWeights.emplace_back();
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
    outEdges[u].swap(neighbours);
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
    assert(n == G.upperNodeIdBound());
    #ifdef NETWORKIT_SANITY_CHECKS
    assert(G.checkConsistency());
    #endif
    
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
    G.setEdgeCount(unsafe, numberOfEdges(G));
    assert(n == G.upperNodeIdBound());
    #ifdef NETWORKIT_SANITY_CHECKS
    assert(G.checkConsistency());
    #endif
    G.shrinkToFit();

    reset();

    return G;
}

void GraphBuilder::toGraphDirectSwap(Graph &G) {

    for (index u = 0; u < outEdges.size(); u++){
        for (index j = 0; j < outEdges[u].size(); j++){
            if(weighted){
                G.addPartialEdge(unsafe, u, outEdges[u][j], outEdgeWeights[u][j]);
            } else {
                G.addPartialEdge(unsafe, u, outEdges[u][j]);
            }
        }
    }

    if(directed){
        for (index u = 0; u < inEdges.size(); u++){
            for (index j = 0; j < inEdges[u].size(); j++){
                if(weighted){
                    G.addPartialInEdge(unsafe, u, inEdges[u][j], inEdgeWeights[u][j]);
                } else {
                    G.addPartialInEdge(unsafe, u, inEdges[u][j]);
                }
            }
        }
    }

    if (!directed) {
        G.setNumberOfSelfLoops(unsafe, selfloops);
    } else if (selfloops % 2 == 0) {
        G.setNumberOfSelfLoops(unsafe, selfloops / 2);
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

    using Adjacencylists = std::vector<std::vector<node>>;
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
                    G.addPartialEdge(unsafe, v, u, outEdgeWeights[v][i]);
                    inWeightsPerThread[tid][u].push_back(outEdgeWeights[v][i]);
                } else {
                    G.addPartialEdge(unsafe, v, u);
                }
            }
            if (u == v) {
                numberOfSelfLoopsPerThread[tid]++;
                if (!directed) {
                    inEdgesPerThread[tid][u].push_back(v);
                    if (weighted)
                        inWeightsPerThread[tid][u].push_back(outEdgeWeights[v][i]);
                }
            }
        }
    });

    // step 2
    parallelForNodes([&](node v) {
        count inDeg = 0;
        count outDeg = G.degreeOut(v);
        for (int tid = 0; tid < maxThreads; tid++) {
            inDeg += inEdgesPerThread[tid][v].size();
        }

        assert(inDeg <= n);
        assert(outDeg <= n);
        // allocate enough memory for all edges/weights
        if (directed) {
            G.preallocateDirectedInEdges(v, inDeg);
        } else {
            G.preallocateUndirected(v, outDeg + inDeg);
        }

        // collect 'second' half of the edges
        if (directed) {
            if (weighted) {
                for (int tid = 0; tid < maxThreads; tid++) {
                    for (index j = 0; j < inEdgesPerThread[tid][v].size(); j++){
                        G.addPartialInEdge(unsafe, v, inEdgesPerThread[tid][v][j], inWeightsPerThread[tid][v][j]);
                    }
                }
            } else {
                for (int tid = 0; tid < maxThreads; tid++) {
                    for (index j = 0; j < inEdgesPerThread[tid][v].size(); j++){
                        G.addPartialInEdge(unsafe, v, inEdgesPerThread[tid][v][j]);
                    }
                }
            }
        } else {
            if (weighted) {
                for (int tid = 0; tid < maxThreads; tid++) {
                    for (index j = 0; j < inEdgesPerThread[tid][v].size(); j++){
                        G.addPartialEdge(unsafe, v, inEdgesPerThread[tid][v][j], inWeightsPerThread[tid][v][j]);
                    }
                }
            } else {
                for (int tid = 0; tid < maxThreads; tid++) {
                    for (index j = 0; j < inEdgesPerThread[tid][v].size(); j++){
                        G.addPartialEdge(unsafe, v, inEdgesPerThread[tid][v][j]);
                    }
                }
            }
        }
    });

    count numSelfLoops = 0;
#pragma omp parallel for reduction(+ : numSelfLoops)
    for (omp_index i = 0; i < static_cast<omp_index>(maxThreads); ++i)
        numSelfLoops += numberOfSelfLoopsPerThread[i];
    G.setNumberOfSelfLoops(unsafe, numSelfLoops);

}

void GraphBuilder::toGraphSequential(Graph &G) {
    G.removeAllEdges();
    
    std::vector<count> missingEdgesCounts(n, 0);
    count numberOfSelfLoops = 0;

    // 'first' half of the edges
    for (index u = 0; u < outEdges.size(); u++){
        for (index j = 0; j < outEdges[u].size(); j++){
            if(weighted) {
                G.addPartialEdge(unsafe, u, outEdges[u][j], outEdgeWeights[u][j]);
            } else {
                G.addPartialEdge(unsafe, u, outEdges[u][j]);
            }
        }
    }
    
    std::vector<count> outDeg(G.upperNodeIdBound());
    parallelForNodes([&](node v) { outDeg[v] = G.degreeOut(v); });

    // count missing edges for each node
    G.forNodes([&](node v) {
        // increase count of incoming edges for all neighbors
        G.forNeighborsOf(v , [&](node u) {
            if (directed || u != v) {
                ++missingEdgesCounts[u];
            }
            if (u == v) {
                // self loops don't need to be added again
                // but we need to count them
                ++numberOfSelfLoops;
            }
        });
    });

    // 'second' half the edges
    if (directed) {
        // directed: outEdges is complete, missing half edges are the inEdges
        // missingEdgesCounts are our inDegrees
        std::vector<count> inDeg = missingEdgesCounts;

        // reserve the exact amount of space needed
        G.forNodes([&](node v) {
            G.preallocateDirectedInEdges(v, inDeg[v]);
        });

        // copy values
        G.forNodes([&](node v) {
            for (index i = 0; i < outDeg[v]; i++) {
                auto nodeWithWeight = G.getIthNeighborWithWeight(unsafe, v, i);
                node u = nodeWithWeight.first;
                if (weighted) {
                    edgeweight ew = nodeWithWeight.second;
                    G.addPartialInEdge(unsafe, u, v, ew);
                } else {
                    G.addPartialInEdge(unsafe, u, v);
                }
            }
        });
    } else {
        // undirected: so far each edge is just saved at one node
        // add it to the other node as well

        // reserve the exact amount of space needed
        G.forNodes([&](node v) {
            G.preallocateUndirected(v, outDeg[v] + missingEdgesCounts[v]);
        });

        // copy values
        G.forNodes([&](node v) {
            // the first outDeg[v] edges in G.outEdges[v] are the first half edges
            // we are adding after outDeg[v]
            for (index i = 0; i < outDeg[v]; i++) {
                auto nodeWithWeight = G.getIthNeighborWithWeight(unsafe, v, i);
                node u = nodeWithWeight.first;
                if (u != v) {
                    if (weighted) {
                        edgeweight ew = nodeWithWeight.second;
                        G.addPartialEdge(unsafe, u, v, ew);
                    } else {
                        G.addPartialEdge(unsafe, u, v);
                    }
                } else {
                    // ignore self loops here
                }
            }
        });
    }

    G.setNumberOfSelfLoops(unsafe, numberOfSelfLoops);
}

count GraphBuilder::numberOfEdges(const Graph &G) {
    count m = 0;
#pragma omp parallel for reduction(+ : m) if (n > (1 << 20))
    for (omp_index v = 0; v < static_cast<omp_index>(G.upperNodeIdBound()); v++) {
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
