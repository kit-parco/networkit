/*
 * BidirectionalBFS.cpp
 *
 *  Created on: 14.06.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#include <limits>

#include <networkit/distance/BidirectionalBFS.hpp>

namespace NetworKit {

void BidirectionalBFS::run() {
    if (G->isWeighted())
        WARN("Treating the graph as unweighted!");

    stDist = 0;
    if (source == target) {
        hasRun = true;
        return;
    }

    init();
    visited.resize(G->upperNodeIdBound(), ts);
    if (ts++ == 128) {
        ts = 1;
        std::fill(visited.begin(), visited.end(), 0);
    }
    if (storePred) {
        pred[source] = source;
        pred[target] = target;
    }

    // Two queues per BFS: sQueue(next) for the exploration from the source node, tQueue(next) for
    // the exploration from the target node.
    sQueue.push(source);
    tQueue.push(target);

    // Mark source as visited by the ball we grow from the source
    visited[source] = ts;
    // Mark the target as visited the ball we grow from the target
    visited[target] = ts + ballMask;

    bool stop = false;
    node s = none, t = none;

    // Expands a ball of a level; idx is the ball index (either 0 or 1<<7)
    auto expand = [&](std::queue<node> &q, std::queue<node> &q_, uint8_t idx) {
        auto visitEdge = [&](node u, node v) {
            if (ts != (visited[v] & ~ballMask)) {
                q_.push(v);
                visited[v] = ts + idx;
                if (storePred) {
                    pred[v] = u;
                }
            } else if ((visited[v] & ballMask) != idx) {
                // Balls met
                s = idx ? v : u;
                t = idx ? u : v;
                stop = true;
            }
        };

        do {
            node u = q.front();
            q.pop();

            if (!G->isDirected() || !idx) {
                // Expanding from source
                for (node v : G->neighborRange(u)) {
                    visitEdge(u, v);
                    if (stop)
                        break;
                }
            } else {
                // Expanding from target
                for (node v : G->inNeighborRange(u)) {
                    visitEdge(u, v);
                    if (stop)
                        break;
                }
            }

        } while (!q.empty());
        std::swap(q, q_);
    };

    do {
        ++stDist;
        if (sQueue.size() <= tQueue.size())
            expand(sQueue, sQueueNext, 0);
        else
            expand(tQueue, tQueueNext, ballMask);
    } while (!stop && !sQueue.empty() && !tQueue.empty());

    // Balls did not meet, source cannot reach target
    if (!stop)
        stDist = std::numeric_limits<count>::max();
    else if (storePred) {
        assert(s != none);
        assert(t != none);

        // Reverse predecessors to target
        node t1 = pred[t];
        pred[t] = s;
        while (t1 != target) {
            node t2 = pred[t1];
            pred[t1] = t;
            t = t1;
            t1 = t2;
        }
        if (t1 != t)
            pred[t1] = t;

        buildPath();
    }

    hasRun = true;
}
} // namespace NetworKit
