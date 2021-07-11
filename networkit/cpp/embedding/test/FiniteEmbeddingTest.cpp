/*
 * FiniteEmbeddingTest.cpp
 *
 *  Created on: 16.10.2020
 *      Author: Klaus Ahrens
 */

#include <iomanip>
#include <iostream>

#include <gtest/gtest.h>

#include <networkit/embedding/Node2Vec.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class FiniteEmbeddingTest : public testing::Test {};

using Embeddings = std::vector<std::vector<float>>;

bool allFinite(const Embeddings &e) {
    for (auto &emb : e) {
        for (auto f : emb) {
            if (std::isinf(f)) {
                return false;
            }
        }
    }
    return true;
}

TEST_F(FiniteEmbeddingTest, testFiniteEmbeddingOnSmallGraph) {

    NetworKit::METISGraphReader reader;

    NetworKit::Graph graph = reader.read("input/karate.graph");

    auto algo = NetworKit::Node2Vec(graph, .3, 1.3);

    algo.run();

    auto features = algo.getFeatures();

    EXPECT_TRUE(allFinite(features));
}

} // namespace NetworKit
