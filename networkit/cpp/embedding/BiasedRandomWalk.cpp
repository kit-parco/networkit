/*
 * BiasedRandomWalk.cpp
 *
 *  Created on: 03.07.2020
 *      Author: Klaus Ahrens  <ahrens@informatik.hu-berlin.de>
 *
 *  adapted and reimplemented from node2vec
 *  part of snap [https://github.com/snap-stanford/snap]
 *  Copyright (c) 2007-2019, Jure Leskovec (under BSD license)
 *
 *  see [https://arxiv.org/pdf/1607.00653v1.pdf]
 *
 */

#include <algorithm>
#include <utility>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/Graph.hpp>

#include "AliasSampler.hpp"
#include "BiasedRandomWalk.hpp"

namespace NetworKit {
namespace Embedding {

using Walk = BiasedRandomWalk::Walk;
using AllWalks = BiasedRandomWalk::AllWalks;

BiasedRandomWalk::BiasedRandomWalk(const Graph *graph)
    : graph(graph), graphData(std::make_unique<GraphData>(graph->numberOfNodes())) {
    auto nn = graph->numberOfNodes();
    index2node.resize(nn);
    // pre-allocate unordered maps:
    graph->forNodes([&](node v) {
        auto degv = graph->degreeOut(v);
        for (auto n : graph->neighborRange(v)) {
            // init index2node:
            index2node[v].push_back(n);
            graphData->data[v][n] = AliasSampler(degv);
        }
    });
}

void BiasedRandomWalk::preprocessNode(node t, double paramP, double paramQ) {

    // for node t
    NeighborSet tNbrs; // Neighbors of t
    for (auto i : graph->neighborRange(t)) {
        tNbrs.insert(i);
    }
    // for each t-neighbor v:
    for (auto v : graph->neighborRange(t)) {
        double pSum = 0;
        std::vector<float> pTable; // Probability distribution table
        // for each v-neighbor
        graph->forNeighborsOf(v, [&](node x, edgeweight weight) {
            if (x == t) {
                pTable.push_back(weight / paramP);
                pSum += weight / paramP;
            } else if (tNbrs.count(x)) {
                pTable.push_back(weight);
                pSum += weight;
            } else {
                pTable.push_back(weight / paramQ);
                pSum += weight / paramQ;
            }
        });
        // Normalizing table
        float pfSum = (float)pSum;
        std::for_each(pTable.begin(), pTable.end(), [pfSum](float &p) { p /= pfSum; });

        graphData->data[v].find(t)->second.unigram(pTable);
    }
}

// Preprocess transition probabilities for each path t->v->x
void BiasedRandomWalk::preprocessTransitionProbs(double paramP, double paramQ) {
    graph->balancedParallelForNodes([&](node t) { preprocessNode(t, paramP, paramQ); });
}

// Simulates a random walk
Walk BiasedRandomWalk::oneWalk(node start, count walkLen) {
    Walk walk(walkLen);
    count nr = 0;
    walk[nr++] = start;
    node src = start;

    if (walkLen == 1) {
        return walk;
    }
    if (graph->degreeOut(start) == 0) {
        walk.resize(1); // shorten walk to 1
        return walk;
    }
    auto nn = graph->degree(start);

    node randomNeighbor = index2node[start][Aux::Random::index(nn)];

    walk[nr++] = randomNeighbor;
    node dst = randomNeighbor;

    while (nr < walkLen) {
        if (graph->degreeOut(dst) == 0) {
            walk.resize(nr); // shorten walk to nr
            return walk;
        }
        NeighborMap &map = graphData->data[dst];
        AliasSampler &as = map[src];
        node next = index2node[dst][as.sample()];
        walk[nr++] = next;
        src = dst;
        dst = next;
    }
    return walk;
}

struct WalkData {
    // assuming contiguous node numbers 0..n-1, they serve as indices
    AllWalks data;
    count walkLength;
    count walksPerNode;

    WalkData(count nn, count wl, count wpn)
        : data{nn * wpn, Walk(wl)}, walkLength(wl), walksPerNode(wpn) {}
};

/// Simulates walks from every node and writes it into Walks vector
AllWalks BiasedRandomWalk::doWalks(count walkLen, count walksPerNode) {
    auto nn = graph->numberOfNodes();

    WalkData walkData(nn, walkLen, walksPerNode);
    std::vector<node> shuffled(graph->nodeRange().begin(), graph->nodeRange().end());

    for (index c = 0; c < walksPerNode; ++c) {
        std::shuffle(shuffled.begin(), shuffled.end(), Aux::Random::getURNG());
        graph->balancedParallelForNodes([&](node i) {
            auto v = shuffled[i];
            Walk thisWalk = oneWalk(v, walkLen);
            walkData.data[c * nn + i] = std::move(thisWalk);
        });
    }

    graphData.reset();

    return walkData.data;
}

} // namespace Embedding
} // namespace NetworKit
