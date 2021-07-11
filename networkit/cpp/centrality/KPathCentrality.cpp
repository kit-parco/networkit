// no-networkit-format
/*
 * KPathCentrality.cpp
 *
 *  Created on: 05.10.2014
 *      Author: nemes
 */

#include <cassert>
#include <stack>

#include <networkit/centrality/KPathCentrality.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

KPathCentrality::KPathCentrality(const Graph& G, double alpha, count k) : Centrality(G, false, false) {
    if (alpha >= -0.5 && alpha <= 0.5) {
        this->alpha = alpha;
    } else {
        throw std::runtime_error("alpha must lie in interval [-0.5, 0.5]");
    }
    if (k == 0) {
        this->k = log(static_cast<double>(G.numberOfNodes() + G.numberOfEdges()));
    } else if (k >0) {
        this->k = k;
    } else {
        throw std::runtime_error("k must be an integer");
    }
}

void KPathCentrality::run() {
    count z = G.upperNodeIdBound();
    double n = G.numberOfNodes();
    scoreData.clear();
    scoreData.resize(z);

    std::vector<count> counter;
    std::vector<bool> explored;

    counter.assign(z, 0);
    explored.assign(z, false);

    count t = 2.0 * k * k * pow(n, 1 - 2 * alpha) * log(n);
    std::stack<node> stack;
    node v = none;

    for (index i = 1; i <= t; i++) { // FIXME: int -> count
        node s = GraphTools::randomNode(G);
        auto l = Aux::Random::integer(1, k);
        explored[s] = true;
        stack.push(s);
        count j = 1;

        while (j <= l) {
            edgeweight sum = 0;
            std::vector<node> neighbours;
            neighbours.clear();
            std::vector<edgeweight> weights;
            weights.clear();
            G.forNeighborsOf(s, [&](node u, edgeweight ew) {
                if (!explored[u]) {
                    neighbours.push_back(u);
                    weights.push_back(1/ew);
                    sum += 1/ew;
                }
            });
            if (neighbours.empty()) {
                break;
            }
            if (G.isWeighted()) {
                double random = Aux::Random::real(0, sum);
                for (index x = 0; x < weights.size(); x++) {
                    if (random < weights[x]) {
                        v = neighbours[x];
                        break;
                    }
                    random -= weights[x];
                }
            } else {
                v = neighbours[Aux::Random::integer(0,neighbours.size() - 1)];
            }
            assert(v != none);
            explored[v] = true;
            stack.push(v);
            counter[v]++;
            s = v;
            j++;
        }
        while (!stack.empty()) {
            v = stack.top();
            stack.pop();
            explored[v] = false;
        }
    }

    G.forNodes([&](node v) {
        scoreData[v] = k * n * (static_cast<double>(counter[v]) / t);
    });

    hasRun = true;
}


} /* namespace NetworKit */
