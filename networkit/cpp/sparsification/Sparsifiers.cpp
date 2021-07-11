// no-networkit-format
/*
 * Sparsifiers.cpp
 *
 *  Created on: 24.07.2014
 *      Author: Gerd Lindner
 */

#include <networkit/sparsification/Sparsifiers.hpp>
#include <networkit/edgescores/TriangleEdgeScore.hpp>
#include <networkit/edgescores/PrefixJaccardScore.hpp>
#include <networkit/sparsification/SimmelianOverlapScore.hpp>
#include <networkit/sparsification/MultiscaleScore.hpp>
#include <networkit/sparsification/LocalSimilarityScore.hpp>
#include <networkit/sparsification/RandomEdgeScore.hpp>
#include <networkit/sparsification/GlobalThresholdFilter.hpp>

#include <networkit/auxiliary/Random.hpp>

namespace NetworKit {

Sparsifier::Sparsifier(const Graph& inputGraph) : inputGraph(inputGraph) {
}

Graph Sparsifier::getGraph() {
    if (!hasOutput) throw std::runtime_error("Error: run must be called first");

    hasOutput = false;
    return std::move(outputGraph);
}


SimmelianSparsifierNonParametric::SimmelianSparsifierNonParametric(const Graph& graph, double threshold) :
        Sparsifier(graph), threshold(threshold) {}

void SimmelianSparsifierNonParametric::run() {
    TriangleEdgeScore triangleEdgeScore(inputGraph);
    triangleEdgeScore.run();
    std::vector<count> triangles = triangleEdgeScore.scores();

    PrefixJaccardScore<count> jaccardScore(inputGraph, triangles);
    jaccardScore.run();
    std::vector<double> jaccard = jaccardScore.scores();

    GlobalThresholdFilter filter(inputGraph, jaccard, threshold, true);
    outputGraph = filter.calculate();
    hasOutput = true;
}


SimmelianSparsifierParametric::SimmelianSparsifierParametric(const Graph& graph, int maxRank, int minOverlap) :
        Sparsifier(graph), maxRank(maxRank), minOverlap(minOverlap) {}

void SimmelianSparsifierParametric::run() {
    TriangleEdgeScore triangleEdgeScore(inputGraph);
    triangleEdgeScore.run();
    std::vector<count> triangles = triangleEdgeScore.scores();

    SimmelianOverlapScore overlapScore(inputGraph, triangles, maxRank);
    overlapScore.run();
    std::vector<double> overlap = overlapScore.scores();

    GlobalThresholdFilter filter(inputGraph, overlap, minOverlap, true);
    outputGraph = filter.calculate();
    hasOutput = true;
}


MultiscaleSparsifier::MultiscaleSparsifier(const Graph& graph, double alpha) :
        Sparsifier(graph), alpha(alpha) {}

void MultiscaleSparsifier::run() {
    std::vector<double> weight(inputGraph.upperEdgeIdBound());
    inputGraph.forEdges([&](node, node, edgeweight w, edgeid eid) {
        weight[eid] = w;
    });

    MultiscaleScore multiscaleScorer(inputGraph, weight);
    multiscaleScorer.run();
    std::vector<double> multiscale = multiscaleScorer.scores();

    GlobalThresholdFilter filter(inputGraph, multiscale, alpha, true);
    outputGraph = filter.calculate();
    hasOutput = true;
}


LocalSimilaritySparsifier::LocalSimilaritySparsifier(const Graph& graph, double e) :
        Sparsifier(graph), e(e) {}

void LocalSimilaritySparsifier::run() {
    TriangleEdgeScore triangleEdgeScore(inputGraph);
    triangleEdgeScore.run();
    std::vector<count> triangles = triangleEdgeScore.scores();

    LocalSimilarityScore localSimScore(inputGraph, triangles);
    localSimScore.run();
    std::vector<double> minExponent = localSimScore.scores();

    GlobalThresholdFilter filter(inputGraph, minExponent, e, true);
    outputGraph = filter.calculate();
    hasOutput = true;
}

SimmelianMultiscaleSparsifier::SimmelianMultiscaleSparsifier(const Graph& graph, double alpha) :
        Sparsifier(graph), alpha(alpha) {}

void SimmelianMultiscaleSparsifier::run() {
    TriangleEdgeScore triangleEdgeScore(inputGraph);
    triangleEdgeScore.run();
    std::vector<count> triangles = triangleEdgeScore.scores();
    std::vector<double> triangles_d = std::vector<double>(triangles.begin(), triangles.end());

    MultiscaleScore multiscaleScorer (inputGraph, triangles_d);
    multiscaleScorer.run();
    std::vector<double> multiscale = multiscaleScorer.scores();

    GlobalThresholdFilter filter(inputGraph, multiscale, alpha, true);
    outputGraph = filter.calculate();
    hasOutput = true;
}

RandomSparsifier::RandomSparsifier(const Graph& graph, double ratio) :
        Sparsifier(graph), ratio(ratio) {}

void RandomSparsifier::run() {
    RandomEdgeScore randomScorer (inputGraph);
    randomScorer.run();
    std::vector<double> random = randomScorer.scores();

    GlobalThresholdFilter filter(inputGraph, random, ratio, true);
    outputGraph = filter.calculate();
    hasOutput = true;
}

} /* namespace NetworKit */
