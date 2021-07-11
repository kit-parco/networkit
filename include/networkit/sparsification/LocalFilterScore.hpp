// no-networkit-format
/*
 * LocalFilterScore.hpp
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#ifndef NETWORKIT_SPARSIFICATION_LOCAL_FILTER_SCORE_HPP_
#define NETWORKIT_SPARSIFICATION_LOCAL_FILTER_SCORE_HPP_

#include <atomic>
#include <memory>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

/**
 * Local filtering edge scoring. Edges with high score are more important.
 *
 * Edges are ranked locally, the top d^e (logarithmic, default) or 1+e*(d-1) edges (non-logarithmic) are kept.
 * For equal attribute values, neighbors of low degree are preferred.
 */
template<typename InType>
class LocalFilterScore final : public EdgeScore<double> {

public:
    /**
     * Initialize the local edge filtering score.
     *
     * @param G The graph for which the score shall be.
     * @param attribute The input attribute according to which the edges shall be filtered locally.
     * @param logarithmic If the score shall be logarithmic in the rank (then d^e edges are kept). Linear otherwise.
     */
    LocalFilterScore(const Graph& G, const std::vector< InType > &attribute, bool logarithmic = true) :
        EdgeScore<double>(G), attribute(&attribute), logarithmic(logarithmic) {}

    /**
     * Execute the algorithm.
     */
    void run() override {
        if (!G->hasEdgeIds()) {
            throw std::runtime_error("edges have not been indexed - call indexEdges first");
        }

        /*
        * For each edge, we calculate the minimum required sparsification exponent e
        * such that the edge is contained in the sparse graph.
        */

        std::unique_ptr<std::atomic<double>[]> sparsificationExp(new std::atomic<double>[G->upperEdgeIdBound()]{});

        G->balancedParallelForNodes([&](node i) {
            count d = G->degree(i);

            /*
             * The top d^e edges (sorted by similarity in descending order)
             * are to be kept in the sparse graph.
             */

            std::vector<edgeid> neighbors;
            neighbors.reserve(d);
            G->forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
                neighbors.emplace_back(eid);
            });

            std::sort(neighbors.begin(), neighbors.end(), [&](const edgeid& e1, const edgeid& e2) {
                return (*attribute)[e1] > (*attribute)[e2];
            });

            count rank = 0;
            count numSame = 1;
            InType oldValue = std::numeric_limits<InType>::lowest();

            for (edgeid eid : neighbors) {
                if ((*attribute)[eid] != oldValue) {
                    rank += numSame;
                    numSame = 1;
                } else {
                    ++numSame;
                }

                double e = 1.0;

                if (d > 1) {
                    if (logarithmic) {
                        e = 1.0 - log(rank) / log(d);
                    } else {
                        e = 1.0 - (rank-1) * 1.0 / (d - 1); // Keep top 1 + e * (d-1) edges
                    }
                }

                Aux::Parallel::atomic_max(sparsificationExp[eid], e);
            }

        });

        scoreData.clear();
        scoreData.resize(G->upperEdgeIdBound());

        #pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(scoreData.size()); ++i) {
            scoreData[i] = sparsificationExp[i];
        }

        hasRun = true;
    }

    double score(node u, node v) override {
        throw std::runtime_error("Not implemented: Use scores() instead.");
    }

    double score(edgeid eid) override {
        throw std::runtime_error("Not implemented: Use scores() instead.");
    }

private:
    const std::vector<InType>* attribute;
    bool logarithmic;

};

} // namespace NetworKit

#endif // NETWORKIT_SPARSIFICATION_LOCAL_FILTER_SCORE_HPP_
