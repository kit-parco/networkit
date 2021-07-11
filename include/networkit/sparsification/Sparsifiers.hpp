// no-networkit-format
/*
 * Sparsifiers.h
 *
 *  Created on: 23.07.2014
 *      Author: Gerd Lindner
 */

#ifndef NETWORKIT_SPARSIFICATION_SPARSIFIERS_HPP_
#define NETWORKIT_SPARSIFICATION_SPARSIFIERS_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * In this file, edge score calculators and edge score filters are combined
 * into different sparsification algorithms.
 */

/**
 * Abstract base class for Sparsifiers.
 */
class Sparsifier {

public:

    Sparsifier(const Graph& inputGraph);

    virtual ~Sparsifier() = default;

    /**
     * REQ: Needs to fill outputGraph and set hasOutput to true.
     */
    virtual void run() = 0;

    Graph getGraph();

protected:
    const Graph& inputGraph;
    Graph outputGraph;
    bool hasOutput;

};


/**
 * Imlementation of the non-parametric variant of Simmelian Backbones,
 * as introduced by Nick et al.
 */
class SimmelianSparsifierNonParametric : public Sparsifier {

public:
    /**
     * Creates a new instance of the non-parametric (jaccard) variant of the Simmelian Backbone calculator.
     * @param graph			the input graph
     * @param threshold		the jaccard index threshold.
     */
    SimmelianSparsifierNonParametric(const Graph& graph, double threshold);

    void run() override;

private:
    double threshold;

};

 /**
  * Imlementation of the parametric variant (Top-k neighborhood overlap)
  * of Simmelian Backbones, as introduced by Nick et al.
  */
class SimmelianSparsifierParametric : public Sparsifier {

public:
    /**
     * Creates a new instance of the parametric variant of the Simmelian Backbone calculator.
     * @param graph			the input graph
     * @param maxRank 		the maximum rank that is considered for overlap calculation
     * @param minOverlap	the minimum overlap of the top-k neighbors for an edge to be
                                contained in the sparsified graph.
     */
    SimmelianSparsifierParametric(const Graph& graph, int maxRank, int minOverlap);

    void run() override;

private:
    int maxRank;
    int minOverlap;

};

/**
 * Implementation of the Multiscale Backbone, as introduced by Serrano et al.
 */
class MultiscaleSparsifier : public Sparsifier {

public:
    /**
     * Creates a new instance of the Multiscale Backbone calculator.
     * @param graph			the input graph
     * @param alpha 		the probability threshold
     */
    MultiscaleSparsifier(const Graph& graph, double alpha);

    void run() override;

private:
    double alpha;

};

/**
 * Local Similarity Sparsification as introduced by Satuluri et al.
 */
class LocalSimilaritySparsifier : public Sparsifier {

public:
    /**
     * Creates a new instance of the Local Similarity sparsifier.
     * @param graph			the input graph
     * @param e				the threshold value
     */
    LocalSimilaritySparsifier(const Graph& graph, double e);

    void run() override;

private:
    double e;

};

/**
 * Multiscale backbone using simmelianness as input weight.
 */
class SimmelianMultiscaleSparsifier : public Sparsifier {

public:
    /**
     * Creates a new instance of the Simmelian Multiscale Backbone calculator
     * @param graph			the input graph
     * @param alpha			the threshold value for multiscale filtering
     */
    SimmelianMultiscaleSparsifier(const Graph& graph, double alpha);

    void run() override;

private:
    double alpha;

};

/**
* Produces sparsified graphs that contain approximately a given percentage
* of edges of the original graph. The edges are selected unformly at random.
*/
class RandomSparsifier : public Sparsifier {

public:
    /**
    * Creates a new instance of the Random Sparsifier.
    * @param graph			the input graph
    * @param ratio			edge ratio in [0,1] to be kept in the sparse graph.
    */
    RandomSparsifier(const Graph& graph, double ratio);

    void run() override;

private:
    double ratio;

};

} /* namespace NetworKit */
#endif // NETWORKIT_SPARSIFICATION_SPARSIFIERS_HPP_
