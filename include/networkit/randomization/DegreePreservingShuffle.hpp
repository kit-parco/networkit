/*
 * DegreePreservingShuffle.hpp
 *
 *  Created on: 21.08.2018
 *      Author: Manuel Penschuck <networkit@manuel.jetzt>
 */

#ifndef NETWORKIT_RANDOMIZATION_DEGREE_PRESERVING_SHUFFLE_HPP_
#define NETWORKIT_RANDOMIZATION_DEGREE_PRESERVING_SHUFFLE_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * 	Implementation of the preprocessing step proposed in
 *  "Smaller Universes for Uniform Sampling of 0,1-matrices with fixed row and column sums"
 *  by Annabell Berger, Corrie Jacobien Carstens [https://arxiv.org/abs/1803.02624]
 *
 *  The algorithms randomizes a graph without changing its topology simply
 *  by renaming nodes. For any degree d (in case of an directed graph it's a degree pair)
 *  consider the set X_d of node ids which have this degree. Then shuffle the ids in X_d.
 *
 *  Hence the algorithm satisfies: For all x in Ginput:
 *   i)  Ginput.degreeIn(x) = Goutput.degreeIn(x)
 *   ii) Ginput.degreeOut(x) = Goutput.degreeOut(x)
 *
 *  The authors argue that applying this preprocessing step before executing (Global)Curveball
 *  leads to a faster mixing time.
 *
 *  @node If you want to use it as a preprocessing step to GlobalCurveball, it's more
 *  efficient to set degreePreservingShufflePreprocessing in GlobalCurveball's constructor.
 */
class DegreePreservingShuffle final : public Algorithm {
public:
    DegreePreservingShuffle() = delete;

    /**
     * Instantiate a DegreePreservingShuffle object
     *
     * @param G Input graph that will be shuffled
     */
    explicit DegreePreservingShuffle(const Graph &G);

    ~DegreePreservingShuffle() override;

    /**
     * Execute trades as configured in the constructor. The algorithm is parallel.
     * @warning This function has to be called exactly once before invoking getGraph()
     */
    void run() final;

    /**
     * Returns a shuffled copy of the input graph.
     * @warning Invoke run() before calling this function.
     */
    Graph getGraph() const;

    /**
     * Returns a reference to the permutation used for shuffling,
     * with permutation[old] = new.
     *
     * @warning Invoke run() before calling this function.
     */
    const std::vector<node> &getPermutation() const noexcept { return permutation; }

    std::string TLX_DEPRECATED(toString() const final);

    bool TLX_DEPRECATED(isParallel() const final) { return true; }

    static Graph shuffleGraph(const Graph &input) {
        DegreePreservingShuffle algo(input);
        algo.run();
        return algo.getGraph();
    }

private:
    const Graph *G;
    std::vector<node> permutation;
};

} // namespace NetworKit

#endif // NETWORKIT_RANDOMIZATION_DEGREE_PRESERVING_SHUFFLE_HPP_
