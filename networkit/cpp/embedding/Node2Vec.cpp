/*
 *  Node2Vec.cpp
 *
 *
 *  Created on: 29.09.2020
 *      Author: Klaus Ahrens
 *              <ahrens@informatik.hu-berlin.de>
 *
 *  adapted and reimplemented from node2vec
 *  part of snap [https://github.com/snap-stanford/snap]
 *  Copyright (c) 2007-2019, Jure Leskovec (under BSD license)
 *
 *  see [https://arxiv.org/pdf/1607.00653v1.pdf]
 *
 */

#include <vector>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/embedding/Node2Vec.hpp>
#include <networkit/graph/Graph.hpp>

#include "AliasSampler.hpp"
#include "BiasedRandomWalk.hpp"
#include "LearnEmbedding.hpp"

namespace NetworKit {

Node2Vec::Node2Vec(const Graph &G, double P, double Q, count L, count N, count D)
    : G(&G), P(P), Q(Q), L(L), N(N), D(D) {
    if (G.numberOfNodes() != G.upperNodeIdBound()) {
        throw std::runtime_error("The node ids of the graph must be continuous.");
    }
}

void Node2Vec::run() {
    using namespace Embedding;

    Aux::SignalHandler handler;

    TRACE("preprocess transition probabilities ...");

    BiasedRandomWalk brw(G);
    brw.preprocessTransitionProbs(P, Q);
    handler.assureRunning();

    TRACE("do biased walks ...");
    auto walks = brw.doWalks(L, N);
    handler.assureRunning();

    TRACE("learn embeddings ...");
    count winSize = 10;
    count iterations = 1;
    features = learnEmbeddings(walks, G->numberOfNodes(), D, winSize, iterations);
    handler.assureRunning();

    hasRun = true;
}

const std::vector<std::vector<float>> &Node2Vec::getFeatures() {
    assureFinished();
    return features;
}

} /* namespace NetworKit */
