/*
 * AliasSampler.hpp
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

#ifndef ALIASSAMPLER_HPP
#define ALIASSAMPLER_HPP

#include <vector>

#include <networkit/Globals.hpp>

namespace NetworKit {
namespace Embedding {

struct AliasSampler {
    std::vector<index> K;
    std::vector<float> U;

public:
    // compute unigram table using alias sampling method
    void unigram(const std::vector<float> &probs);

    AliasSampler(count degree = 0)
        : K(degree), U(degree) {} // allowing as mapped_type in unordered_maps
    index sample();
};

} // namespace Embedding
} // namespace NetworKit

#endif // ALIASSAMPLER_HPP
