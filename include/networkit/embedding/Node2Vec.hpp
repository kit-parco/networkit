/*
 *  Node2Vec.hpp
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
 */

#ifndef NETWORKIT_EMBEDDING_NODE2_VEC_HPP_
#define NETWORKIT_EMBEDDING_NODE2_VEC_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup embedding
 *
 * Feature extraction by node2vec algorithm
 */
class Node2Vec final : public Algorithm {

public:
    /**
     * Creates the Node2Vec class for @a G
     *
     * @param G     The graph.
     * @param P     Walk Return parameter (stay local).
     * @param Q     Walk In-Out parameter (drift away).
     * @param L     Walk length.
     * @param N     Walk count.
     * @param D     Dimension of embedding vectors.
     */
    Node2Vec(const Graph &G, double P = 1, double Q = 1, count L = 80, count N = 10, count D = 128);

    ~Node2Vec() override = default;

    /**
     This method computes all node embeddings.
    */
    void run() override;

    /**
     This method returns a vector that contains feature vectors for all nodes
    */
    const std::vector<std::vector<float>> &getFeatures();

private:
    // The graph
    const Graph *G;
    // Walk parameter P
    double P;
    // Walk parameter Q
    double Q;
    // Walk length L
    count L;
    // Walk count N
    count N;
    // Feature Dimension D
    count D;

    std::vector<std::vector<float>> features;
};

} /* namespace NetworKit */

#endif // NETWORKIT_EMBEDDING_NODE2_VEC_HPP_
