// no-networkit-format
/*
 * CoreDecomposition.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Lukas Barth, David Wei, Christian Staudt
 *  Inplace change on Jun 26, 2015 by Henning Meyerhenke
 */

#include <omp.h>
#include <set>

#include <networkit/centrality/DegreeCentrality.hpp>
#include <networkit/centrality/CoreDecomposition.hpp>

namespace NetworKit {

CoreDecomposition::CoreDecomposition(const Graph& G, bool normalized, bool enforceBucketQueueAlgorithm, bool storeNodeOrder) :
        Centrality(G, normalized), maxCore(0), enforceBucketQueueAlgorithm(enforceBucketQueueAlgorithm),
        storeNodeOrder(storeNodeOrder)
{
    if (G.numberOfSelfLoops()) throw std::runtime_error("Core Decomposition implementation does not support graphs with self-loops. Call Graph.removeSelfLoops() first.");
    if (storeNodeOrder) this->enforceBucketQueueAlgorithm = true;
    canRunInParallel = (! enforceBucketQueueAlgorithm && (G.numberOfNodes() == G.upperNodeIdBound()));
}

void CoreDecomposition::run() {
    if (G.isDirected() || enforceBucketQueueAlgorithm) {
        runWithBucketQueues();
    }
    else {
        runWithParK();
    }

    if (normalized) {
        DegreeCentrality deg(G);
        deg.run();
        auto degrees = deg.scores();
        count maxDeg = *std::max_element(degrees.begin(), degrees.end());
        G.parallelForNodes([&](node u) {
            scoreData[u] = scoreData[u] / maxDeg;
        });
    }
}

void CoreDecomposition::runWithParK() {
    count z = G.upperNodeIdBound();
    scoreData.resize(z); // TODO: move to base class

    count nUnprocessed = G.numberOfNodes();
    std::vector<node> curr; // currently processed nodes
    std::vector<node> next; // nodes to be processed next
    std::vector<char> active(z,0);
    index level = 0; // current level
    count size = 0;  // number of nodes currently processed

    // fill in degrees
    std::vector<count> degrees(z);
    G.parallelForNodes([&](node u) {
        degrees[u] = G.degree(u);
        active[u] = 1;
    });

    // main loop
    while (nUnprocessed > 0) {
        // find nodes with degree == current level
        if (! canRunInParallel || z <= 256) {
            scan(level, degrees, curr);
        }
        else {
            scanParallel(level, degrees, curr, active);
        }

        // process such nodes in curr
        size = curr.size();
        while (size > 0) {
            nUnprocessed -= size;
#ifndef NETWORKIT_OMP2
            if (! canRunInParallel || size <= 256) {
                processSublevel(level, degrees, curr, next);
            }
            else {
                processSublevelParallel(level, degrees, curr, next, active);
            }
#else
            processSublevel(level, degrees, curr, next);
#endif
            std::swap(curr, next);
            size = curr.size();
            next.clear();
        }
        ++level;
    }

    maxCore = level-1;
    hasRun = true;
}

void CoreDecomposition::scan(index level, const std::vector<count>& degrees,
        std::vector<node>& curr)
{
    G.forNodes([&](node u) {
        if (degrees[u] == level) {
            curr.push_back(u);
        }
    });
}

void CoreDecomposition::scanParallel(index level, const std::vector<count>& degrees,
        std::vector<node>& curr, std::vector<char>& active)
{
    const count z = G.upperNodeIdBound();
    std::vector<std::vector<node>> next(omp_get_max_threads());
    curr.clear();

#pragma omp parallel for schedule(guided)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        if (active[u] && degrees[u] == level) {
            auto tid = omp_get_thread_num();
            next[tid].push_back(u);
        }
    }
    for (auto& n : next) {
        curr.insert(curr.end(),n.begin(),n.end());
    }
}

void CoreDecomposition::processSublevel(index level,
        std::vector<count>& degrees, const std::vector<node>& curr,
        std::vector<node>& next)
{
    // check for each neighbor of vertices in curr if their updated degree reaches level;
    // if so, process them next
    for (auto u: curr) {
        scoreData[u] = level;
        G.forNeighborsOf(u, [&](node v) {
            if (degrees[v] > level) {
                degrees[v]--;
                if (degrees[v] == level) {
                    next.push_back(v);
                }
            }
        });
    }
}

#ifndef NETWORKIT_OMP2
void CoreDecomposition::processSublevelParallel(index level,
        std::vector<count>& degrees, const std::vector<node>& curr,
        std::vector<node>& next, std::vector<char>& active)
{
    // check for each neighbor of vertices in curr if their updated degree reaches level;
    // if so, process them next

    const count size = curr.size();
    std::vector<std::vector<node>> localNext(omp_get_max_threads());

#pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(size); ++i) {
        node u = curr[i];
        active[u] = 0;
        scoreData[u] = level;
        G.forNeighborsOf(u, [&](node v) {
            if (degrees[v] > level) {
                index tmp;

#pragma omp atomic capture
                tmp = --degrees[v];

                // ensure that neighbor is inserted exactly once if necessary
                if (tmp == level) { // error in external publication on ParK fixed here
                    auto tid = omp_get_thread_num();
                    localNext[tid].push_back(v);
                }
            }
        });
    }
    for (auto& n : localNext) {
        next.insert(next.end(), n.begin(), n.end());
    }
}
#endif // NETWORKIT_OMP2

void CoreDecomposition::runWithBucketQueues() {
    /* Main data structure: buckets of nodes indexed by their remaining degree. */
    index z = G.upperNodeIdBound();
    std::vector<node> queue(G.numberOfNodes());
    std::vector<index> nodePtr(z);
    std::vector<index> degreeBegin(G.numberOfNodes());
    std::vector<count> degree(z);       // tracks degree during algo

    bool directed = G.isDirected();

    /* Bucket sort  by degree */
    /* 1) bucket sizes */
    G.forNodes([&](node u) {
        count deg = G.degree(u);

        if (directed) {
            deg += G.degreeIn(u);
        }

        degree[u] = deg;
        ++degreeBegin[deg];
    });

    index sum = 0; // 2) exclusive in-place prefix sum
    for (index i = 0; i < degreeBegin.size(); ++i) {
        index tmp = degreeBegin[i];
        degreeBegin[i] = sum;
        sum += tmp;
    }

    /* 3) Sort nodes/place them in queue */
    G.forNodes([&](node u) {
        count deg = degree[u];
        index pos = degreeBegin[deg];
        ++degreeBegin[deg];
        queue[pos] = u;
        nodePtr[u] = pos;
    });

    /* 4) restore exclusive prefix sum */
    index tmp = 0; // move all one forward
    for (index i = 0; i < degreeBegin.size(); ++i) {
        std::swap(tmp, degreeBegin[i]);
    }

    /* Current core and and computed scoreData values. */
    index core = 0;
    scoreData.clear();
    scoreData.resize(z);

    /* Main loop: Successively "remove" nodes by setting them not alive after processing them. */
    for (index i = 0; i < G.numberOfNodes(); ++i) {
        node u = queue[i];
        core = std::max(core, degree[u]); // core is maximum of all previously seen degrees

        scoreData[u] = core;

        /* Remove a neighbor by decreasing its degree and changing its position in the queue */
        auto removeNeighbor = [&](node v) {
            if (nodePtr[v] > i) { // only nodes that are after the current node need to be considered
                // adjust the degree
                count oldDeg = degree[v];
                --degree[v];
                count newDeg = oldDeg - 1;

                // Degrees smaller than the current degree can be before the current position
                // Correct those that we need. Note that as we decrease degrees only by one
                // all degreeBegin values larger than oldDeg will have been adjusted already.
                if (degreeBegin[oldDeg] <= i) {
                    degreeBegin[oldDeg] = i + 1;
                }

                if (degreeBegin[newDeg] <= i) {
                    degreeBegin[newDeg] = i + 1;
                }

                // Swap v with beginning of the current bucket.
                index oldPos = nodePtr[v];
                index newPos = degreeBegin[oldDeg];
                node nodeToSwap = queue[newPos];
                std::swap(queue[oldPos], queue[newPos]);
                std::swap(nodePtr[nodeToSwap], nodePtr[v]);

                // Move bucket border, v is now in the previous bucket, i.e. the bucket of its new degree
                ++degreeBegin[oldDeg];
            }
        };

        /* Remove u and its incident edges. */
        G.forNeighborsOf(u, removeNeighbor);
        if (directed) {
            /* graph is directed */
            G.forInNeighborsOf(u, removeNeighbor);
        }
    }

    maxCore = core;

    if (storeNodeOrder) {
        nodeOrder = std::move(queue);
    }

    hasRun = true;
}


Cover CoreDecomposition::getCover() const {
    if (! hasRun) throw std::runtime_error("call run method first");
    // initialize Cover
    index z = G.upperNodeIdBound();
    Cover coverData;
    if (coverData.numberOfElements() != z) {
        coverData = Cover(z);
        coverData.setUpperBound(z);
    }
    // enter values from scoreData into coverData
    G.forNodes([&](node u) {
        index k = 0;
        while (scoreData[u] >= k) {
            coverData.addToSubset((index) k, (index) u);
            ++k;
        }
    });
    return coverData;
}

Partition CoreDecomposition::getPartition() const {
    if (! hasRun) throw std::runtime_error("call run method first");
    // initialize Partition
    index z = G.upperNodeIdBound();
    Partition shellData;
    if (shellData.numberOfElements() != z) {
        shellData = Partition(z);
        shellData.allToSingletons();
    }
    // enter values from scoreData into shellData
    G.forNodes([&](node u){
        shellData.moveToSubset((index) scoreData[u], (index) u);
    });
    return shellData;
}

const std::vector<node>& CoreDecomposition::getNodeOrder() const {
    if (!storeNodeOrder) throw std::runtime_error("The node order was not stored. Make sure you set storeNodeOrder to true");
    assureFinished();

    return nodeOrder;
}

index CoreDecomposition::maxCoreNumber() const {
    if (! hasRun) throw std::runtime_error("call run method first");
    return maxCore;
}

double CoreDecomposition::maximum() {
    return static_cast<double>(G.numberOfNodes() - 1);
}

} /* namespace NetworKit */
