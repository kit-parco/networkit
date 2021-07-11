/*
 * LearnEmbeddings.cpp
 *
 *  Created on: 06.07.2020
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
#include <vector>

#include <networkit/auxiliary/Random.hpp>

#include "AliasSampler.hpp"
#include "BiasedRandomWalk.hpp"
#include "LearnEmbedding.hpp"

// Code from https://github.com/nicholas-leonard/word2vec/blob/master/word2vec.c
// Customized for SNAP and node2vec

/*
The kernel algorithm in trainModel contains several threading issues:

as in the original code there are concurrent reads and writes on
alpha, wordCntAll, synPos and synNeg !

Nevertheless, missing some writes should not harm the functionality of the
algorithm, but add a further (indefinite) source of randomness.
*/

namespace NetworKit {
namespace Embedding {

// Number of negative samples. Value taken from original word2vec code.
constexpr int negSamN = 5;

constexpr int maxExp = 6;

// Learning parameters for SGD. Values taken from original word2vec code.
constexpr double startAlpha = 0.025;
constexpr int refreshAlphaCount = 10000;
constexpr double minAlpha = startAlpha * 0.0001;
constexpr double nsPwr = 0.75;
/*
The exponent used to shape the negative sampling distribution.
A value of 1.0 samples exactly in proportion to the frequencies,
0.0 samples all words equally, while a negative value samples
low-frequency words more than high-frequency words.
The popular default value of 0.75 was chosen by the original Word2Vec paper.
*/

using Vocab = std::vector<count>;

Vocab learnVocab(const AllWalks &allWalks, count nn) {
    Vocab vocab(nn);
    for (auto &walk : allWalks) { // for every walk starting from this node
        for (auto &node : walk) { // for every node in this walk
            ++vocab[node];
        }
    }
    return vocab;
}

// Precompute unigram table using alias sampling method
AliasSampler vocabSampler(const Vocab &vocab) {
    double trainWordsPow = 0;
    auto sz = vocab.size();
    std::vector<float> probV(sz); // NOT ... {sz}
    for (count i = 0; i < sz; ++i) {
        double e = std::pow(vocab[i], nsPwr);
        probV[i] = e;
        trainWordsPow += e;
    }
    std::transform(probV.begin(), probV.end(), probV.begin(),
                   [&trainWordsPow](float f) { return f / trainWordsPow; });

    AliasSampler as(sz);
    as.unigram(probV);
    return as;
}

// Initialize negative embeddings
Embeddings initNegEmb(const Vocab &vocab, int dimensions) {
    return Embeddings{vocab.size(), std::vector<float>(dimensions)};
}

// Initialize positive embeddings
Embeddings initPosEmb(const Vocab &vocab, int dimensions) {
    Embeddings emb{vocab.size(), std::vector<float>(dimensions)};
#ifndef NETWORKIT_OMP2
#pragma omp parallel for schedule(dynamic)
#endif // NETWORKIT_OMP2
    for (omp_index i = 0; i < static_cast<omp_index>(emb.size()); ++i) {
        std::generate(emb[i].begin(), emb[i].end(),
                      [&]() { return (Aux::Random::real() - 0.5) / dimensions; });
    }
    return emb;
}

void trainModel(ModelData &model, count walkNr) {
    AllWalks &allWalks = model.allWalks;
    count &wordCntAll = model.wordCntAll;
    AliasSampler &as = model.as;
    count allWords = model.allWords;
    count iterations = model.iterations;
    count dimensions = model.dimensions;
    count winSize = model.winSize;
    double &alpha = model.alpha;
    Embeddings &synPos = model.synPos;
    Embeddings &synNeg = model.synNeg;

    using Walk = BiasedRandomWalk::Walk;

    Walk &thisWalk = allWalks[walkNr];

    for (index wordI = 0; wordI < thisWalk.size(); ++wordI) {
        count localWordCntAll;
#ifndef NETWORKIT_OMP2
#pragma omp atomic read
#endif // NETWORKIT_OMP2
        localWordCntAll = wordCntAll;
        if (localWordCntAll % refreshAlphaCount == 0) {
            double newAlpha =
                startAlpha * (1 - localWordCntAll / ((double)(iterations * allWords) + 1.0));
#ifndef NETWORKIT_OMP2
#pragma omp atomic write
#endif // NETWORKIT_OMP2
            alpha = (newAlpha < minAlpha) ? minAlpha : newAlpha;
        }

        node word = thisWalk[wordI];
        std::vector<float> eV(dimensions);
        index offset = Aux::Random::index(winSize);

        for (index a = offset; a < winSize * 2 + 1 - offset; ++a) {
            if (a == winSize) {
                continue;
            }
            if (wordI + a < winSize) {
                continue;
            }
            if (wordI + a >= winSize + thisWalk.size()) {
                continue;
            }
            count currWordI = wordI + a - winSize;
            count currWord = thisWalk[currWordI];
            // negative sampling
            for (index j = 0; j < negSamN + 1; ++j) {
                node target;
                int label;
                if (j == 0) {
                    target = word;
                    label = 1;
                } else {
                    target = as.sample();
                    if (target == word) {
                        continue;
                    }
                    label = 0;
                }
                double product = 0;
                for (index i = 0; i < dimensions; ++i) {
                    product += synPos[currWord][i] * synNeg[target][i];
                }
                double grad; // Gradient multiplied by learning rate
                if (product > maxExp) {
                    grad = (label - 1) * alpha;
                } else if (product < -maxExp) {
                    grad = label * alpha;
                } else {
                    double exp = std::exp(product);
                    grad = (label - 1 + 1 / (1 + exp)) * alpha;
                }

                for (index i = 0; i < dimensions; ++i) {
                    float sNti;
#ifndef NETWORKIT_OMP2
#pragma omp atomic read
#endif // NETWORKIT_OMP2
                    sNti = synNeg[target][i];
                    eV[i] += grad * sNti;
#ifndef NETWORKIT_OMP2
#pragma omp atomic write
#endif // NETWORKIT_OMP2
                    synNeg[target][i] = sNti + grad * synPos[currWord][i];
                }
            }

            for (index i = 0; i < dimensions; ++i) {
                float sPci;
#ifndef NETWORKIT_OMP2
#pragma omp atomic read
#endif // NETWORKIT_OMP2
                sPci = synPos[currWord][i];

#ifndef NETWORKIT_OMP2
#pragma omp atomic write
#endif // NETWORKIT_OMP2
                synPos[currWord][i] = sPci + eV[i];
            }
        }
#ifndef NETWORKIT_OMP2
#pragma omp atomic write
#endif // NETWORKIT_OMP2
        wordCntAll = localWordCntAll + 1;
    }
}

Embeddings learnEmbeddings(AllWalks &allWalks, count nn, count dimensions, count winSize,
                           count iterations) {
    count walks = allWalks.size();           // number of walks
    count nodesPerWalk = allWalks[0].size(); // number of nodes per walk

    std::vector<node> renumbered(nn, none);
    std::vector<node> original(nn, none);
    count nodes = 0;
    for (index i = 0; i < walks; ++i) {
        for (index j = 0; j < nodesPerWalk; ++j) {
            node &nd = allWalks[i][j];
            if (renumbered[nd] != none) {
                nd = renumbered[nd];
            } else {
                renumbered[nd] = nodes;
                original[nodes] = nd;
                nd = nodes++;
            }
        }
    }

    std::vector<count> vocab = learnVocab(allWalks, nodes);

    Embeddings embeddings(nodes, std::vector<float>(dimensions));

    Embeddings synNeg = initNegEmb(vocab, dimensions);
    Embeddings synPos = initPosEmb(vocab, dimensions);
    AliasSampler as = vocabSampler(vocab);

    double alpha = startAlpha; // learning rate
    count wordCntAll = 0;

    ModelData model(allWalks, dimensions, winSize, iterations, as, wordCntAll, alpha, synNeg,
                    synPos);

    for (index iterCnt = 0; iterCnt < iterations; ++iterCnt) {
#ifndef NETWORKIT_OMP2
#pragma omp parallel for schedule(dynamic)
#endif // NETWORKIT_OMP2
        for (omp_index run = 0; run < static_cast<omp_index>(walks); ++run) {
            trainModel(model, run);
        }
    }

    for (index i = 0; i < nodes; ++i) {
        embeddings[original[i]] = synPos[i];
    }

    allWalks.clear();
    return embeddings;
}

} // namespace Embedding
} // namespace NetworKit
