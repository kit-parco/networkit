/*
 * AliasSampler.cpp
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

#include <random>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Random.hpp>

#include "AliasSampler.hpp"

namespace NetworKit {
namespace Embedding {

void AliasSampler::unigram(const std::vector<float> &probs) {
    auto N = probs.size();

    std::vector<index> underV;
    std::vector<index> overV;
    for (index i = 0; i < N; ++i) {
        U[i] = probs[i] * N;
        if (U[i] < 1) {
            underV.push_back(i);
        } else {
            overV.push_back(i);
        }
    }
    while (!underV.empty() && !overV.empty()) {
        auto small = underV.back();
        auto large = overV.back();
        underV.pop_back();
        overV.pop_back();
        K[small] = large;
        U[large] = (U[large] + U[small]) - 1;
        if (U[large] < 1) {
            underV.push_back(large);
        } else {
            overV.push_back(large);
        }
    }
    for (index i : underV)
        U[i] = 1;
    for (index i : overV)
        U[i] = 1;
}

index AliasSampler::sample() {
    index x = static_cast<index>(Aux::Random::index(K.size()));
    return Aux::Random::real() < U[x] ? x : K[x];
}

} // namespace Embedding
} // namespace NetworKit
