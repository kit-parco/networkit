/*
 * EdgeSwitching.hpp
 *
 *  Created on: 19.10.2019
 *      Author:  Manuel Penschuck <networkit@manuel.jetzt>
 */

#ifndef NETWORKIT_RANDOMIZATION_EDGE_SWITCHING_HPP_
#define NETWORKIT_RANDOMIZATION_EDGE_SWITCHING_HPP_

#include <random>
#include <string>
#include <utility>

#include <networkit/Globals.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * EdgeSwitching and EdgeSwitchingInPlace implement the same algorithm -- while the former copies
 * the input graph provided and is more versatile, the in-place implementation only keeps a
 * reference to the input graph.
 *
 * The Edge Switching Markov Chain ["The markov chain simulation method for generating connected
 * power law random graphs", Mihail and Zegura] perturbates simple directed or undirected graphs
 * while preserving their degrees. In each step, two edges are selected uniformly at random, and
 * their endpoints exchanged. Swaps that introduce multi-edges or self-loops are rejected WITHOUT
 * replacement -- this is necessary to allow uniform sampling [see "Switching edges to randomize
 * networks: what goes wrong and how to fix it", Carstens and Horadam]. The number of successful
 * swaps can be queried using getNumberOfAffectedEdges()/2.
 *
 * @note In general, simple edge switching does not yield a uniform distribution on simple DIRECTED
 * graphs because the orientation of directed triangles cannot be changed. Using
 * DegreePreservingShuffle as a preprocessing step overcomes this limitation. The preprocessing can
 * also jump-start the perturbation process, yielding to potentially faster mixing.
 *
 * It's recommended to implement the preprocessing step in the calling code or to use copying
 * EdgeSwitching implementation.
 */
class EdgeSwitchingInPlace : public Algorithm {
public:
    /**
     * Constructs an EdgeSwitch algorithm that will change the input graph IN-PLACE.
     * @warning For directed graphs preprocessing with DegreePreservingShuffle is necessary. Either
     * do it manually, or use another constructor.
     */
    EdgeSwitchingInPlace(Graph &G, double numberOfSwitchesPerEdge = 10.0) : graph(G) {
        setNumberOfSwitchesPerEdge(numberOfSwitchesPerEdge);
    }

    ~EdgeSwitchingInPlace() override = default;

    /**
     * Attempts to carry out numberOfSwitchesPerEdge * numberOfEdges() edge swaps (i.e., rejected
     * swaps are counted as well, see Algorithm's description why)
     * @note This function can be called several times. The graph's topology may not be changed
     * externally between two invocations.
     */
    void run() override;

    /// Return a reference to the perturbed graph
    Graph &getGraph() const { return graph; }

    /**
     * Returns twice the total number of non-rejected swaps carried out during calls to run().
     * It hence has the same semantics as Curveball::getNumberOfAffectedEdges().
     * @note If run() is called multiple time, the sum over all invocations is returned.
     */
    count getNumberOfAffectedEdges() const noexcept { return 2 * numberOfSwapsPerformed; }

    /// Return (average) number of switches per edge that will be executed on next run
    double getNumberOfSwitchesPerEdge() const noexcept { return numberOfSwitchesPerEdge; }

    /// Modify (average) number of switches per edge that will be executed on next run
    void setNumberOfSwitchesPerEdge(double x);

protected:
    Graph &graph;
    std::discrete_distribution<node>
        degreeDistribution; ///< Return node x with probability proportional to degree(x)

    double numberOfSwitchesPerEdge;
    count numberOfSwapsPerformed{0};
};

/**
 * EdgeSwitching and EdgeSwitchingInPlace implement the same algorithm -- while the former copies
 * the input graph provided and is more versatile, the in-place implementation only keeps a
 * reference to the input graph.
 *
 * The Edge Switching Markov Chain ["The markov chain simulation method for generating connected
 * power law random graphs", Mihail and Zegura] perturbates simple directed or undirected graphs
 * while preserving their degrees. In each step, two edges are selected uniformly at random, and
 * their endpoints exchanged. Swaps that introduce multi-edges or self-loops are rejected WITHOUT
 * replacement -- this is necessary to allow uniform sampling [see "Switching edges to randomize
 * networks: what goes wrong and how to fix it", Carstens and Horadam]. The number of successful
 * swaps can be queried using getNumberOfAffectedEdges()/2.
 *
 * @note In general, simple edge switching does not yield a uniform distribution on simple DIRECTED
 * graphs because the orientation of directed triangles cannot be changed. Using
 * DegreePreservingShuffle as a preprocessing step overcomes this limitation. The preprocessing can
 * also jump-start the perturbation process, yielding to potentially faster mixing.
 *
 * It's recommended to implement the preprocessing step in the calling code or to use copying
 * EdgeSwitching implementation.
 */
class EdgeSwitching : public Algorithm {
public:
    /// Constructs an EdgeSwitch algorithm that contains a COPY of the input graph.
    explicit EdgeSwitching(const Graph &G, double numberOfSwitchesPerEdge = 10.0,
                           bool degreePreservingShufflePreprocessing = true);

    ~EdgeSwitching() override = default;

    /**
     * Attempts to carry out numberOfSwitchesPerEdge * numberOfEdges() edge swaps (i.e., rejected
     * swaps are counted as well, see Algorithm's description why)
     */
    void run() override { inPlaceAlgorithm.run(); }

    /// Return a reference to the perturbed graph
    const Graph &getGraph() const { return ownedGraph; }

    /**
     * Move graph owned by the algorithm out.
     * @warning Do not call run() after calling moveGraph()
     */
    Graph moveGraph() { return std::move(ownedGraph); }

    /**
     * Returns twice the total number of non-rejected swaps carried out during calls to run().
     * It hence has the same semantics as Curveball::getNumberOfAffectedEdges().
     * @note If run() is called multiple time, the sum over all invocations is returned.
     */
    count getNumberOfAffectedEdges() const noexcept {
        return inPlaceAlgorithm.getNumberOfAffectedEdges();
    }

    /// Return (average) number of switches per edge that will be executed on next run
    double getNumberOfSwitchesPerEdge() const noexcept {
        return inPlaceAlgorithm.getNumberOfSwitchesPerEdge();
    }

    /// Modify (average) number of switches per edge that will be executed on next run
    void setNumberOfSwitchesPerEdge(double x) { inPlaceAlgorithm.setNumberOfSwitchesPerEdge(x); }

private:
    Graph ownedGraph;
    EdgeSwitchingInPlace inPlaceAlgorithm;
};

} // namespace NetworKit

#endif // NETWORKIT_RANDOMIZATION_EDGE_SWITCHING_HPP_
