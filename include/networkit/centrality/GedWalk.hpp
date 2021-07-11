/* GedWalk.hpp
 *
 *  Created on: 18.6.2018
 *     Authors: Eugenio Angriman <angrimae@hu-berlin.de>
 *              Alexander van der Grinten <avdgrinten@hu-berlin.de>
 */

#ifndef NETWORKIT_CENTRALITY_GED_WALK_HPP_
#define NETWORKIT_CENTRALITY_GED_WALK_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

namespace NetworKit {

class GedWalk final : public Algorithm {

    using walks = double;

public:
    enum GreedyStrategy { lazy, stochastic };

    enum BoundStrategy { no, spectral, geometric, adaptiveGeometric };

    double stocEpsilon = 0.1;

private:
    const Graph *G;
    const count k;
    const double epsilon;
    double alpha;
    const BoundStrategy boundStrategy;
    GreedyStrategy greedyStrategy;
    double spectralDelta;
    double degOutMax;
    double degInMax;
    std::vector<double> alphas;

    // Size of pathsIn / pathsOut and related vectors.
    count allocatedLevels = 15;

    double sigmaMax;

    // Similar to gainScore, gainBound and gainW, but for the entire group.
    // These values are always maintained exactly.
    double groupScore; // Score of the group.
    walks groupW;      // Number of walks on the last level.
    walks groupBound;  // Upper bound on the score of the group.

    walks graphW;

    count nLevels = 2;
    std::vector<node> group;
    std::vector<unsigned char> inGroup;
    std::vector<unsigned char> isExact;
    std::vector<unsigned char> nodesToPick;
    std::vector<std::vector<walks>> pathsIn, pathsOut;
    std::vector<std::vector<walks>> pathsHit, pathsMiss;

    // Unless isExact[x], these values might be upper bounds on the true values.
    std::vector<double> gainScore; // Marginal score of the group.
    std::vector<walks> gainW;      // Number of marginal walks on the last level.
    std::vector<double> gainBound; // Upper bound on the marginal score.

    struct CompareScore {
    public:
        CompareScore(const std::vector<double> &score) : score(&score) {}
        bool operator()(node x, node y) const { return (*score)[x] > (*score)[y]; }

    private:
        const std::vector<double> *score;
    };

    struct EvaluationResult {
        double score;
        walks w;
    };

    tlx::d_ary_addressable_int_heap<node, 2, CompareScore> scoreQ{CompareScore(gainScore)};
    tlx::d_ary_addressable_int_heap<node, 2, CompareScore> boundQ{CompareScore(gainBound)};

    void init();
    double computeGamma();
    void estimateGains();
    void computeMarginalGain(node z);
    void maximizeGain();
    bool separateNodes();
    double computeSigmaMax() const;
    void fillPQs();
    void updateAlphas();
    EvaluationResult evaluateGraph();
    EvaluationResult evaluateGroup();

public:
    /**
     * Finds a group of `k` vertices with at least ((1 - 1/e) * opt - epsilon) GedWalk centrality
     * score, where opt is the highest possible score. The algorithm is based on the paper "Group
     * Centrality Maximization for Large-scale Graphs", Angriman et al., ALENEX20. It implements two
     * independent greedy strategies (lazy and stochastic). Furthermore, it allows to compute the
     * GedWalk score of a given set of nodes.
     *
     * @param G A (weakly) connected graph.
     * @param k The desired group size.
     * @param epsilon Precision of the algorithm.
     * @param alpha Exponent to compute the GedWalk score.
     * @param bs Bound strategy to compute the GedWalk bounds, default: BoundStrategy::geometric.
     * @param gs Greedy strategy to be used (lazy or stochastic), default: GreedyStrategy::lazy.
     * @param spectralDelta Delta to be used for the spectral bound.
     *
     * Note: if the graph is weighted, all weights should be in (0, 1].
     */
    GedWalk(const Graph &G, count k = 1, double initEpsilon = 0.1, double alpha = -1.,
            BoundStrategy bs = BoundStrategy::geometric, GreedyStrategy gs = GreedyStrategy::lazy,
            double spectralDelta = 0.5);

    /**
     * Run the algorithm.
     */
    void run() override;

    /**
     * Return the approximate GedWalk score of the computed group.
     *
     * @return The approximate score of the computed group.
     */
    double getApproximateScore() const {
        assureFinished();
        return groupScore;
    }

    /**
     * Return the GedWalk score of the input group.
     *
     * @param first/last The range that contains the input group.
     * @param scoreEpsilon The precision of the score to be computed.
     */
    template <class InputIt>
    double scoreOfGroup(InputIt first, InputIt last, double scoreEpsilon = .1);

    /**
     * Return the computed group.
     *
     * @return The computed group.
     */
    std::vector<node> groupMaxGedWalk() const {
        assureFinished();
        return group;
    }

    const std::vector<double> &boundOfMarginalGains() const { return gainBound; }
};

template <class InputIt>
double GedWalk::scoreOfGroup(InputIt first, InputIt last, double scoreEpsilon) {
    if (boundStrategy == BoundStrategy::spectral) {
        assert(alpha * static_cast<double>(sigmaMax) < 1.);
    } else if (boundStrategy == BoundStrategy::geometric) {
        assert(alpha * static_cast<double>(degInMax) < 1.);
    } else {
        assert(boundStrategy == BoundStrategy::adaptiveGeometric);
        assert(alpha * static_cast<double>(degOutMax + degInMax) < 1.);
    }

    std::fill(pathsHit[0].begin(), pathsHit[0].end(), walks{0});
    std::fill(pathsMiss[0].begin(), pathsMiss[0].end(), walks{1});
    std::fill(inGroup.begin(), inGroup.end(), static_cast<unsigned char>(0));
    while (first != last) {
        const auto u = *first;
        inGroup[u] = 1;
        pathsHit[0][u] = 1.0;
        pathsMiss[0][u] = 0.0;
        ++first;
    }

    graphW = evaluateGraph().w;

    do {
        const auto result = evaluateGroup();
        groupScore = result.score;
        groupW = result.w;
        if (boundStrategy == BoundStrategy::spectral) {
            const double gamma = sqrt(G->numberOfNodes()) * (sigmaMax / (1 - alpha * sigmaMax));
            groupBound = result.score + alphas[nLevels + 1] * gamma * graphW;
        } else if (boundStrategy == BoundStrategy::geometric) {
            const double gamma = (degInMax / (1 - alpha * degInMax));
            groupBound = result.score + alphas[nLevels + 1] * gamma * graphW;
        } else {
            assert(boundStrategy == BoundStrategy::adaptiveGeometric);
            groupBound = result.score + alphas[nLevels + 1] * computeGamma() * result.w;
        }

        DEBUG("score is ", groupScore, ", bound is ", groupBound);
        if (groupBound < groupScore + scoreEpsilon)
            return groupScore;

        DEBUG("increasing path length to ", nLevels + 1);
        if (nLevels + 1 >= allocatedLevels) {
            DEBUG("allocating ", nLevels + 2, " GedWalk levels");
            while (allocatedLevels < nLevels + 2) {
                pathsIn.emplace_back((G->upperNodeIdBound()));
                pathsOut.emplace_back((G->upperNodeIdBound()));
                pathsHit.emplace_back((G->upperNodeIdBound()));
                pathsMiss.emplace_back((G->upperNodeIdBound()));
                ++allocatedLevels;
            }

            updateAlphas();
        }

        ++nLevels;
        graphW = evaluateGraph().w;
    } while (true);
}

} // namespace NetworKit
#endif // NETWORKIT_CENTRALITY_GED_WALK_HPP_
