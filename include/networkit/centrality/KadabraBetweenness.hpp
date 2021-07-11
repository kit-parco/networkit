// no-networkit-format
/*
 * KadabraBetweenness.hpp
 *
 * Created on: 18.07.2018
 *    Authors: Eugenio Angriman <angrimae@hu-berlin.de>
 *             Alexander van der Grinten <avdgrinten@hu-berlin.de>
 */

#ifndef NETWORKIT_CENTRALITY_KADABRA_BETWEENNESS_HPP_
#define NETWORKIT_CENTRALITY_KADABRA_BETWEENNESS_HPP_

#include <atomic>
#include <memory>
#include <random>
#include <vector>

#include <networkit/auxiliary/SortedList.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class StateFrame {
  public:
    StateFrame(const count size) { apx.assign(size, 0); };
    count nPairs = 0;
    count epoch = 0;
    std::vector<count> apx;

    void reset(count newEpoch) {
        std::fill(apx.begin(), apx.end(), 0);
        nPairs = 0;
        epoch = newEpoch;
    }
};

class Status {
  public:
      Status(count k);
      const count k;
      std::vector<node> top;
      std::vector<double> approxTop;
      std::vector<bool> finished;
      std::vector<double> bet;
      std::vector<double> errL;
      std::vector<double> errU;
};

class SpSampler {
  private:
    const Graph &G;
    const ConnectedComponents &cc;

  public:
    SpSampler(const Graph &G, const ConnectedComponents &cc);
    void randomPath(StateFrame *curFrame);
    StateFrame *frame;
    std::mt19937_64 rng;
    std::uniform_int_distribution<node> distr;

  private:
    std::vector<uint8_t> timestamp;
    uint8_t globalTS = 1;
    static constexpr uint8_t stampMask = 0x7F;
    static constexpr uint8_t ballMask = 0x80;
    std::vector<count> dist;
    std::vector<count> nPaths;
    std::vector<node> q;
    std::vector<std::pair<node, node>> spEdges;

    inline node randomNode() const;
    void backtrackPath(node source, node target, node start);
    void resetSampler(count endQ);
    count getDegree(const Graph &graph, node y, bool useDegreeIn);
};

/**
 * @ingroup centrality
 */
class KadabraBetweenness : public Algorithm {
  public:
    // See EUROPAR'19 paper for the selection of these parameters.
    unsigned int baseItersPerStep = 1000;
    double itersPerStepExp = 1.33;

    /**
     * Approximation of the betweenness centrality and computation of the top-k
     * nodes with highest betweenness centrality according to the algorithm
     * described in Borassi M., and Natale E. (2016): KADABRA is an ADaptive
     * Algorithm for Betweenness via Random Approximation.
     * Parallel implementation by Van der Grinten A., Angriman E., and
     * Meyerhenke H. (2019): Parallel Adaptive Sampling with almost no
     * Synchronization. ArXiv pre-print (https://arxiv.org/abs/1903.09422).
     *
     * If k = 0 the algorithm approximates the betweenness centrality of all
     * vertices of the graph so that the scores are within an additive error @a
     * err with probability at least (1 - @a delta). Otherwise, the algorithm
     * computes the exact ranking of the top-k nodes with highest betweenness
     * centrality. The algorithm relies on an adaptive random sampling technique
     * of shortest paths and the number of samples in the worst case is w =
     * ((log(D - 2) + log(2/delta))/err^2 samples, where D is the diameter of
     * the graph. Thus, the worst-case performance is O(w * (|E| + |V|)), but
     * performs better in practice.
     *
     * @param G The graph
     * @param err Maximum additive error guaranteed when approximating the
     * betweenness centrality of all nodes.
     * @param delta probability that the values of the betweenness centrality
     * are within the error guarantee.
     * @param k the number of top-k nodes to be computed. Set it to zero to
     * approximate the betweenness centrality of all the nodes.
     * @param unionSample, startFactor algorithm parameters that are
     * automatically chosen.
     */
    KadabraBetweenness(const Graph &G, double err = 0.01, double delta = 0.1,
                       bool deterministic = false, count k = 0, count unionSample = 0,
                       count startFactor = 100);

    /**
     * Executes the Kadabra algorithm.
     */
    void run() override;

    /**
     * @return The ranking of the nodes according to their approximated
     * betweenness centrality.
     */
    std::vector<std::pair<node, double>> ranking() const;

    /**
     * @return Nodes of the graph sorted by their approximated betweenness
     * centrality.
     */
    std::vector<node> topkNodesList() const {
        assureFinished();
        return topkNodes;
    }

    /**
     * @return Sorted list of approximated betweenness centrality scores.
     */
    std::vector<double> topkScoresList() const {
        assureFinished();
        return topkScores;
    }

    /**
     * @return Approximated betweenness centrality score of all the nodes of the
     * graph.
     */
    std::vector<double> scores() const {
        assureFinished();
        return approxSum;
    }

    /**
     * @return Total number of samples.
     */
    count getNumberOfIterations() const {
        assureFinished();
        return nPairs;
    }

    /**
     * @return Upper bound to the number of samples.
     */
    double getOmega() const {
        assureFinished();
        return omega;
    }

    count maxAllocatedFrames() const {
        assureFinished();
        count maxNFrames = 0;
        for (auto nFrames : maxFrames) {
            maxNFrames = std::max(nFrames, maxNFrames);
        }
        return maxNFrames;
    }

  protected:
    const Graph &G;
    const double delta, err;
    const bool deterministic;
    const count k, startFactor;
    count unionSample;
    count nPairs;
    const bool absolute;
    double deltaLMinGuess, deltaUMinGuess, omega;
    std::atomic<int32_t> epochToRead;
    int32_t epochRead;
    count seed0, seed1;
    std::vector<count> maxFrames;

    std::vector<node> topkNodes;
    std::vector<double> topkScores;
    std::vector<std::pair<node, double>> rankingVector;
    std::vector<std::atomic<StateFrame *>> epochFinished;
    std::vector<SpSampler> samplerVec;
    std::unique_ptr<Aux::SortedList> top;
    std::unique_ptr<ConnectedComponents> cc;

    std::vector<double> approxSum;
    std::vector<double> deltaLGuess;
    std::vector<double> deltaUGuess;

    std::atomic<bool> stop;

    void init();
    void computeDeltaGuess();
    void computeBetErr(Status *status, std::vector<double> &bet,
                       std::vector<double> &errL,
                       std::vector<double> &errU) const;
    bool computeFinished(Status *status) const;
    void getStatus(Status *status, bool parallel = false) const;
    void computeApproxParallel(const std::vector<StateFrame> &firstFrames);
    double computeF(double btilde, count iterNum, double deltaL) const;
    double computeG(double btilde, count iterNum, double deltaU) const;
    void fillResult();
    void checkConvergence(Status &status);

    void fillPQ() {
        for (count i = 0; i < G.upperNodeIdBound(); ++i) {
            top->insert(i, approxSum[i]);
        }
    }
};

inline std::vector<std::pair<node, double>>
KadabraBetweenness::ranking() const {
    assureFinished();
    std::vector<std::pair<node, double>> result(topkNodes.size());
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(result.size()); ++i) {
        result[i] = std::make_pair(topkNodes[i], topkScores[i]);
    }
    return result;
}
} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_KADABRA_BETWEENNESS_HPP_
