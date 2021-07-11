/* GedWalk.cpp
 *
 *	Created on: 18.6.2018
 *    Authors: Eugenio Angriman <angrimae@hu-berlin.de>
 *             Alexander van der Grinten <avdgrinten@hu-berlin.de>
 */

#include <cmath>
#include <omp.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/centrality/GedWalk.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <tlx/unused.hpp>

namespace NetworKit {

GedWalk::GedWalk(const Graph &graph, count k, double initEpsilon, double alpha, BoundStrategy bs,
                 GreedyStrategy gs, double spectralDelta)
    : G(&graph), k(k), epsilon(initEpsilon / k), alpha(alpha), boundStrategy(bs),
      greedyStrategy(gs), spectralDelta(spectralDelta) {

    if (!k || k >= G->upperNodeIdBound())
        throw std::runtime_error("Error: k should be between 1 and n-1.");

    if (spectralDelta < 0 || spectralDelta > 1)
        throw std::runtime_error("Error: spectralDelta should be between 0 and 1.");

    if (G->isWeighted())
        G->parallelForEdges([&](node, node, const edgeweight ew) {
            tlx::unused(ew);
            assert(ew >= 0.0 && ew <= 1.0);
        });

    init();
}

double GedWalk::computeSigmaMax() const {
    const count n = G->upperNodeIdBound();

    std::vector<double> vector(n, 1.0);
    std::vector<double> tVector(n, 1.0);
    std::vector<double> oldVector(n, 1.0);

    double value = 0.0, oldValue = 0.0;

    auto converged(
        [&](const double val, const double other) -> bool { return std::abs(val - other) < 0.1; });

    do {
        oldValue = value;

        // Iterate matrix-vector product.
        G->parallelForNodes([&](const node u) {
            tVector[u] = 0.0;
            G->forInEdgesOf(
                u, [&](const node v, const edgeweight ew) { tVector[u] += ew * oldVector[v]; });
        });

        G->parallelForNodes([&](const node u) {
            vector[u] = 0.0;
            G->forEdgesOf(u,
                          [&](const node v, const edgeweight ew) { vector[u] += ew * tVector[v]; });
        });

        // Normalize vector.
        value = 0.0;
        value = G->parallelSumForNodes([&](const node u) { return (vector[u] * vector[u]); });
        value = std::sqrt(value);

        G->parallelForNodes([&](const node u) { vector[u] /= value; });

        std::swap(oldVector, vector);
    } while (!converged(value, oldValue));

    DEBUG("sigmaMax is ", sigmaMax, ", degInMax is ", degInMax);
    return std::sqrt(value);
}

void GedWalk::init() {
    const auto n = G->upperNodeIdBound();

    hasRun = false;
    groupScore = 0.0;
    groupW = 0.0;
    groupBound = 0.0;
    degOutMax = static_cast<double>(GraphTools::maxDegree(*G));
    degInMax = static_cast<double>(GraphTools::maxInDegree(*G));

    inGroup.resize(n, 0);
    isExact.resize(n, 0);

    group.clear();
    group.reserve(k);
    gainScore.resize(n);
    gainBound.resize(n);
    gainW.resize(n);

    pathsIn.resize(allocatedLevels, std::vector<walks>(n));
    pathsOut.resize(allocatedLevels, std::vector<walks>(n));
    pathsHit.resize(allocatedLevels, std::vector<walks>(n));
    pathsMiss.resize(allocatedLevels, std::vector<walks>(n));

    scoreQ.reserve(n);
    boundQ.reserve(n);

    if (greedyStrategy == GreedyStrategy::stochastic) {
        nodesToPick.resize(n);
    }

    // For spectral bound: Compute largest singular value.
    if (boundStrategy == BoundStrategy::spectral) {
        sigmaMax = computeSigmaMax();
    }

    if (alpha <= 0) {
        if (boundStrategy == BoundStrategy::spectral) {
            alpha = spectralDelta / sigmaMax;
        } else if (boundStrategy == BoundStrategy::geometric) {
            alpha = 1.0 / (1.0 + degInMax);
        } else {
            assert(boundStrategy == BoundStrategy::adaptiveGeometric);
            alpha = 1.0 / (1.0 + degOutMax + degInMax);
        }
    }

    updateAlphas();
}

void GedWalk::updateAlphas() {
    // nLevels grows up to allocatedLevels - 1; hence, we need alpha[allocatedLevels]
    // when calculating the bound for nLevels = allocatedLevels - 1.0
    const auto alphasSize = alphas.size();
    alphas.resize(allocatedLevels + 1);
    for (auto i = alphasSize; i < alphas.size(); ++i) {
        alphas[i] = std::pow(alpha, i);
    }
}

double GedWalk::computeGamma() {
    auto degSum = degOutMax + degInMax;
    return (degSum / (1.0 - alpha * degSum));
}

/* Computes an upper bound to the GED-walk score of every node.
 * Only used once to initialize the priority queue.
 */
void GedWalk::estimateGains() {
    // Each node has a 0-path.
    // For nodes in the group, we claim that there is no 0-path. This ensures that we do
    // not consider paths that cross the group in this method.
    std::fill(pathsOut[0].begin(), pathsOut[0].end(), walks{1});

    for (const auto u : group) {
        pathsOut[0][u] = 0.0;
    }

    if (G->isDirected()) {
        std::fill(pathsIn[0].begin(), pathsIn[0].end(), walks{1});
        for (const auto u : group) {
            pathsIn[0][u] = 0.0;
        }
    }

    // Compute the number of in-paths and out-paths at all levels.
    const auto n = G->upperNodeIdBound();
    for (count i = 1; i <= nLevels; ++i) {
#pragma omp parallel for
        for (omp_index u = 0; u < static_cast<omp_index>(n); ++u) {
            walks outSum = 0.0;
            if (!inGroup[u])
                G->forNeighborsOf(u, [&](const node v, const edgeweight ew) {
                    outSum += ew * pathsOut[i - 1][v];
                });
            pathsOut[i][u] = outSum;

            if (G->isDirected()) {
                walks inSum = 0.0;
                if (!inGroup[u])
                    G->forInEdgesOf(u, [&](const node v, const edgeweight ew) {
                        inSum += ew * pathsIn[i - 1][v];
                    });
                pathsIn[i][u] = inSum;
            }
        }
    }

    // Accumulate the number of paths that cross a given nodes.
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
        const auto u = static_cast<node>(i);
        // m is the suffix of the path, while l - m is the prefix.
        // Walks ending in u will have a suffix length = 0.

        // Do the calculation for all levels except for the last one.
        double firstLevelsScore = 0.0;
        for (count l = 1; l < nLevels; ++l) {
            walks nWalks = 0.0;
            for (count m = 0; m <= l; ++m) {
                nWalks += (G->isDirected() ? pathsIn : pathsOut)[l - m][u] * pathsOut[m][u];
            }

            firstLevelsScore += alphas[l] * nWalks;
        }

        // Do the calculation for the last level.
        walks w = 0.0;
        for (count m = 0; m <= nLevels; ++m) {
            w += (G->isDirected() ? pathsIn : pathsOut)[nLevels - m][u] * pathsOut[m][u];
        }

        assert(!inGroup[u] || (!firstLevelsScore && !w));

        const double score = firstLevelsScore + alphas[nLevels] * w;
        double bound;
        if (boundStrategy == BoundStrategy::spectral) {
            const double gamma =
                std::sqrt(G->numberOfNodes()) * (sigmaMax / (1 - alpha * sigmaMax));
            bound = firstLevelsScore + alphas[nLevels] * w + alphas[nLevels + 1] * gamma * graphW;
        } else if (boundStrategy == BoundStrategy::geometric) {
            const double gamma = (degInMax / (1 - alpha * degInMax));
            bound = firstLevelsScore + alphas[nLevels] * w + alphas[nLevels + 1] * gamma * graphW;
        } else {
            assert(boundStrategy == BoundStrategy::adaptiveGeometric);
            bound =
                firstLevelsScore + alphas[nLevels] * w + alphas[nLevels + 1] * computeGamma() * w;
        }

        if (score < gainScore[u]) {
            gainScore[u] = score;
            DEBUG("estimated score for ", u, " is ", score);
        }
        if (bound < gainBound[u]) {
            gainBound[u] = bound;
            DEBUG("estimated bound for ", u, " is ", bound, " (", w, " paths)");
        }
        if (w < gainW[u]) {
            gainW[u] = w;
        }
    }

    scoreQ.update_all();
    boundQ.update_all();
}

// Evaluate GED walk for the group consisting of the entire vertex set.
auto GedWalk::evaluateGraph() -> EvaluationResult {
    // Special case the first level: pathsHit[0] is always 1.
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(G->upperNodeIdBound()); ++i) {
        const auto u = static_cast<node>(i);
        walks pHit = 0.0;
        G->forInEdgesOf(u, [&](node, const edgeweight ew) { pHit += ew; });
        pathsHit[1][u] = pHit;
    }

    // Compute the numbers of paths.
    for (count i = 2; i <= nLevels; ++i) {
#pragma omp parallel for
        for (omp_index j = 0; j < static_cast<omp_index>(G->upperNodeIdBound()); ++j) {
            const auto u = static_cast<node>(j);
            walks pHit = 0;
            G->forInEdgesOf(
                u, [&](const node v, const edgeweight ew) { pHit += ew * pathsHit[i - 1][v]; });
            pathsHit[i][u] = pHit;
        }
    }

    // Aggregate the group score.
    double score = 0.0;
    walks w = 0.0;

#ifndef NETWORKIT_OMP2
#pragma omp parallel for reduction(+ : score, w)
    for (omp_index i = 0; i < static_cast<omp_index>(G->upperNodeIdBound()); ++i) {
        const auto u = static_cast<node>(i);
        double contrib = 0.0;
        for (count j = 1; j <= nLevels; ++j) {
            contrib += pathsHit[j][u] * alphas[j];
        }
        score += contrib;
        w += pathsHit[nLevels][u];
    }
#else
    G->forNodes([&](const node u) {
        double contrib = 0.0;
        for (count i = 1; i <= nLevels; ++i) {
            contrib += pathsHit[i][u] * alphas[i];
        }
        score += contrib;
        w += pathsHit[nLevels][u];
    });
#endif

    return EvaluationResult{score, w};
}

auto GedWalk::evaluateGroup() -> EvaluationResult {
    // Compute the numbers of paths.
    for (count i = 1; i <= nLevels; ++i) {
#pragma omp parallel for
        for (omp_index j = 0; j < static_cast<omp_index>(G->upperNodeIdBound()); ++j) {
            const auto u = static_cast<node>(j);
            walks pHit = 0.0, pMiss = 0.0;
            if (inGroup[u])
                G->forInEdgesOf(u, [&](const node v, const edgeweight ew) {
                    pHit += ew * (pathsMiss[i - 1][v] + pathsHit[i - 1][v]);
                });
            else
                G->forInEdgesOf(u, [&](const node v, const edgeweight ew) {
                    pMiss += ew * pathsMiss[i - 1][v];
                    pHit += ew * pathsHit[i - 1][v];
                });

            pathsHit[i][u] = pHit;
            pathsMiss[i][u] = pMiss;
        }
    }

    // Aggregate the group score.
    double score = 0.0;
    walks w = 0.0;

#ifndef NETWORKIT_OMP2
#pragma omp parallel for reduction(+ : score, w)
    for (omp_index i = 0; i < static_cast<omp_index>(G->upperNodeIdBound()); ++i) {
        const auto u = static_cast<node>(i);
        double contrib = 0.0;
        for (count j = 1; j <= nLevels; ++j) {
            contrib += pathsHit[j][u] * alphas[j];
        }
        score += contrib;
        w += pathsHit[nLevels][u];
    }
#else
    G->forNodes([&](const node u) {
        double contrib = 0.0;
        for (count i = 1; i <= nLevels; ++i) {
            contrib += pathsHit[i][u] * alphas[i];
        }
        score += contrib;
        w += pathsHit[nLevels][u];
    });
#endif

    return EvaluationResult{score, w};
}

void GedWalk::computeMarginalGain(node z) {
    assert(!inGroup[z]);

    // For simplicity
    inGroup[z] = 1;
    pathsHit[0][z] = 1;
    pathsMiss[0][z] = 0;

    const auto result = evaluateGroup();

    // Restoring coherence
    inGroup[z] = 0;
    pathsHit[0][z] = 0;
    pathsMiss[0][z] = 1;

    // Score and bound increase monotonically with the group (monotonicity of GED).
    const auto newGroupScore = result.score;
    const auto newGroupW = result.w;

    assert(newGroupScore >= groupScore - epsilon);

    // New scores must never exceed old bounds.
    const auto newGainScore = newGroupScore - groupScore;
    const auto newGainW = newGroupW - groupW;
    double newGainBound;
    if (boundStrategy == BoundStrategy::geometric) {
        const double gamma = (degInMax / (1 - alpha * degInMax));
        newGainBound = newGainScore + alphas[nLevels + 1] * gamma * graphW;
    } else if (boundStrategy == BoundStrategy::spectral) {
        const double gamma = std::sqrt(G->numberOfNodes()) * (sigmaMax / (1 - alpha * sigmaMax));
        newGainBound = newGainScore + alphas[nLevels + 1] * gamma * graphW;
    } else {
        assert(boundStrategy == BoundStrategy::adaptiveGeometric);
        newGainBound = newGainScore + alphas[nLevels + 1] * computeGamma() * newGainW;
    }

    assert(gainBound[z] >= newGainScore - epsilon);

    // Marginal gain and bounds decrease monotonically with the group (submodularity of GED).
    // (If they don't, we cannot keep the priority queues when updating the group.)
    assert(gainScore[z] >= newGainScore - epsilon);
    assert(gainBound[z] >= newGainBound - epsilon);

    gainScore[z] = newGainScore;
    gainW[z] = newGainW;
    gainBound[z] = newGainBound;

    DEBUG(z, scoreQ.top());

    assert(scoreQ.contains(scoreQ.top()));
    assert(boundQ.contains(z));
    scoreQ.update(z);
    boundQ.update(z);

    isExact[z] = 1;
}

void GedWalk::maximizeGain() {
    assert(!scoreQ.empty());

    // TODO: We can stop earlier if gainBound for an isExact[] vertex is higher than
    //       the gainScore of the top vertex.

    // Find the vertex with highest gainScore.
    count gainIters = 0;
    while (!isExact[scoreQ.top()]) {
        computeMarginalGain(scoreQ.top());
        ++gainIters;
    }
    DEBUG(gainIters, " iterations to find maximal gain");
}

bool GedWalk::separateNodes() {
    auto isSeparated = [&](const node z, const node s) {
        return gainScore[z] >= gainBound[s] - epsilon;
    };

    assert(!scoreQ.empty());
    const auto z = scoreQ.top();
    assert(isExact[z]);

    // To compute bounds, temporarily take out our candidate vertex.
    boundQ.remove(z);
    assert(!boundQ.empty());

    // Find the vertex with highest gainBound.
    node s;
    count separationIters = 0;
    do {
        s = boundQ.top();
        assert(s != z);

        if (isExact[s]) {
            break;
        }

        computeMarginalGain(s);
        ++separationIters;
    } while (isSeparated(z, s));

    DEBUG(separationIters, " iterations for separation");

    // Restore the bound queue.
    boundQ.push(z);

    DEBUG("picking ", z, " yields gain ", gainScore[z], ", and bound for ", s, " is ",
          gainBound[s]);

    return isSeparated(z, s);
}

void GedWalk::fillPQs() {
    if (greedyStrategy == GreedyStrategy::lazy) {
        G->forNodes([&](const node u) {
            scoreQ.update(u);
            boundQ.update(u);
        });
    } else if (greedyStrategy == GreedyStrategy::stochastic) {
        count nSamples =
            std::max(1.0, std::log(1.0 / stocEpsilon) * static_cast<double>(G->upperNodeIdBound())
                              / static_cast<double>(k));
        DEBUG("Picking ", nSamples, " new candidates for stochastic greedy");

        if (nSamples >= G->upperNodeIdBound() - group.size()) {
            WARN("Number of samples is too high, reverting to lazy greedy.");
            greedyStrategy = GreedyStrategy::lazy;
            fillPQs();
            return;
        }

        bool pick; // whether by default nodes should be picked (true), or removed (false)

        // If we want to extract more than n/2 samples, it is faster to select all nodes in the
        // graph as samples, and then remove n - 'nSamples' random nodes.
        auto c = nSamples;
        if (nSamples > (G->upperNodeIdBound() - group.size()) / 2) {
            pick = false;
            c = G->upperNodeIdBound() - nSamples - group.size();
        } else {
            // Otherwise, we iteratively select 'nSamples' random samples
            pick = true;
        }

        DEBUG("n = ", G->upperNodeIdBound(), " |S| = ", group.size(), ", c = ", c);
        DEBUG("Pick = ", pick);
        // Mark all nodes as 'not picked'
        std::fill(nodesToPick.begin(), nodesToPick.end(), static_cast<unsigned char>(!pick));

        for (const auto u : group) {
            nodesToPick[u] = 0;
        }

        scoreQ.clear();
        boundQ.clear();

        auto &gen = Aux::Random::getURNG();
        std::uniform_int_distribution<node> nodesDistr{0, G->upperNodeIdBound()};
        do {
            const auto u = nodesDistr(gen);
            if (G->hasNode(u) && !inGroup[u] && nodesToPick[u] != pick) {
                // Mark node u as 'pickded'
                nodesToPick[u] = pick;
                --c;
            }
        } while (c);

        // Insert all the nodes to be picked in the heaps
        for (count i = 0; i < G->upperNodeIdBound(); ++i) {
            if (nodesToPick[i]) {
                assert(!inGroup[i]);
                scoreQ.push(i);
                boundQ.push(i);
                if (scoreQ.size() == nSamples) {
                    break;
                }
            }
        }

        assert(scoreQ.size() == nSamples);
        assert(boundQ.size() == nSamples);
    }
}

void GedWalk::run() {
    if (boundStrategy == BoundStrategy::spectral) {
        assert(alpha < 1.0 / static_cast<double>(sigmaMax));
    } else if (boundStrategy == BoundStrategy::geometric) {
        assert(alpha < 1.0 / degInMax);
    } else {
        assert(boundStrategy == BoundStrategy::adaptiveGeometric);
        assert(alpha < 1.0 / (degOutMax + degInMax));
    }

    // Should we reset nLevels to zero?
    // Doing this will not change the (qualitative) result of the algorithm but it would let it
    // do a "fresh start" if it already ran before.

    std::fill(inGroup.begin(), inGroup.end(), static_cast<unsigned char>(0));
    std::fill(gainScore.begin(), gainScore.end(), std::numeric_limits<double>::max());
    std::fill(gainBound.begin(), gainBound.end(), std::numeric_limits<double>::max());
    std::fill(pathsHit[0].begin(), pathsHit[0].end(), walks{0});
    std::fill(pathsMiss[0].begin(), pathsMiss[0].end(), walks{1});

    fillPQs();
    graphW = evaluateGraph().w;

    do {
        assert(group.size() < k);
        assert(nLevels < allocatedLevels);

        estimateGains();

        // In each iteration, the inner loop tries to add a vertex to the group.
        // It breaks to the outer loop if nLevels needs to increase.
        do {
            if (group.size() == k) {
                hasRun = true;
                return;
            }

            DEBUG("Group size = ", group.size(), " new iteration");

            // This has to be reset whenever the group changes or nLevles is increased.
            // For simplicity, just reset it here.
            std::fill(isExact.begin(), isExact.end(), static_cast<unsigned char>(0));

            maximizeGain();
            DEBUG("Maximize gain success");
            if (!separateNodes()) {
                break;
            }

            // Add z to the group.
            const auto z = scoreQ.extract_top();
            boundQ.remove(z);

            inGroup[z] = 1;
            group.push_back(z);

            pathsHit[0][z] = 1;
            pathsMiss[0][z] = 0;

            groupScore += gainScore[z];
            groupW += gainW[z];

            if (boundStrategy == BoundStrategy::spectral) {
                const double gamma =
                    std::sqrt(G->numberOfNodes()) * (sigmaMax / (1 - alpha * sigmaMax));
                groupBound = groupScore + alphas[nLevels + 1] * gamma * graphW;
            } else if (boundStrategy == BoundStrategy::geometric) {
                const double gamma = (degInMax / (1 - alpha * degInMax));
                groupBound = groupScore + alphas[nLevels + 1] * gamma * graphW;
            } else {
                assert(boundStrategy == BoundStrategy::adaptiveGeometric);
                groupBound = groupScore + alphas[nLevels + 1] * computeGamma() * groupW;
            }

            if (greedyStrategy == GreedyStrategy::stochastic && group.size() < k) {
                fillPQs();
            }

            DEBUG("adding ", z, " to the group, new score: ", groupScore);
        } while (true);

        DEBUG("increasing path length to ", nLevels + 1);
        if (nLevels + 1 >= allocatedLevels) {
            INFO("allocating ", nLevels + 2, " GedWalk levels");
            while (allocatedLevels < nLevels + 2) {
                pathsIn.emplace_back((G->upperNodeIdBound()));
                pathsOut.emplace_back((G->upperNodeIdBound()));
                pathsHit.emplace_back((G->upperNodeIdBound()));
                pathsMiss.emplace_back((G->upperNodeIdBound()));
                ++allocatedLevels;
            }

            updateAlphas();
        }

#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(G->upperNodeIdBound()); ++i) {
            const auto u = static_cast<node>(i);
            if (inGroup[u])
                continue;

            gainScore[u] = std::numeric_limits<double>::max();
            gainW[u] = std::numeric_limits<double>::max();
        }

        scoreQ.update_all();

        ++nLevels;
        graphW = evaluateGraph().w;

        const auto result = evaluateGroup();
        groupScore = result.score;
        groupW = result.w;

        if (boundStrategy == BoundStrategy::spectral) {
            const double gamma =
                std::sqrt(G->numberOfNodes()) * (sigmaMax / (1 - alpha * sigmaMax));
            groupBound = result.score + alphas[nLevels + 1] * gamma * graphW;
        } else if (boundStrategy == BoundStrategy::geometric) {
            const double gamma = (degInMax / (1 - alpha * degInMax));
            groupBound = result.score + alphas[nLevels + 1] * gamma * graphW;
        } else {
            assert(boundStrategy == BoundStrategy::adaptiveGeometric);
            groupBound = result.score + alphas[nLevels + 1] * computeGamma() * result.w;
        }
    } while (true);
}

} // namespace NetworKit
