// no-networkit-format
/*
 * Luby.cpp
 *
 *  Created on: 27.02.2013
 *      Author: Christian Staudt
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/independentset/Luby.hpp>

namespace NetworKit {

std::vector<bool> Luby::run(const Graph& G) {

    std::vector<bool> I(G.numberOfNodes(), false); // independent set $I = \emptyset$
    std::vector<bool> V(G.numberOfNodes(), true); // instead of pruning the graph, store here whether a node in G is still in G'

    // test if there are no active nodes left (G' is empty)
    auto empty = [&](){
        for (bool a : V) {
            if (a) return false;
        }
        return true;
    };

    // weighted degree filtered for active nodes
    auto weightedDegree = [&](node u){
        double wDeg = 0.0;
        G.forNeighborsOf(u, [&](node v, edgeweight w) {
            if (V[v]) {
                wDeg += w;
            }
        });
        return wDeg;
    };

    auto nodeProbability = [&](node v){
        return 1.0 / (2.0 * weightedDegree(v));
    };



    int i = 0;
    while (! empty()) {
        i += 1;
        DEBUG("Luby iteration #" , i);
        // choose set S - weighted choice of active nodes with probability $1 / 2 \omega(v)$
        std::vector<bool> S(G.numberOfNodes(), false);
        G.parallelForNodes([&](node u){
            if (V[u]) {
                if (Aux::Random::probability() < nodeProbability(u)) {
                    S[u] = true;  // add node to S
                }
            }
        });
        // remove non-independent nodes from S to get S'
        G.parallelForEdges([&](node u, node v) {
            if (u != v) { // exclude self-loops
                if (S[u] & S[v]) { // u and v are not independent (note: S is subset of V')
                    // remove node with smaller degree
                    edgeweight wu = weightedDegree(u);
                    edgeweight wv = weightedDegree(v);
                    if (wu > wv) {
                        S[v] = false;
                    } else if (wv > wu) {
                        S[u] = false;
                    } else { // tie
                        S[v] = false; // arbitrary decision
                    }
                }
            }

        });

        // add S' to I
        G.parallelForNodes([&](node u){
            if (S[u]) {
                I[u] = true;
            }
        });

        // remove S' and all neighboring nodes from V'
        G.parallelForNodes([&](node u){
            if (S[u]) {
                V[u] = false;
                G.forNeighborsOf(u, [&](node v){
                    V[v] = false;
                });
            }
        });

    }

    return I;
}

std::string Luby::toString() const {
    return "Luby";
}

} /* namespace NetworKit */
