// no-networkit-format
/*
 * SparsificationBenchmark.cpp
 *
 *  Created on: 31.07.2014
 *      Author: Gerd Lindner
 */

#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/edgescores/ChibaNishizekiTriangleEdgeScore.hpp>
#include <networkit/edgescores/TriangleEdgeScore.hpp>
#include <networkit/edgescores/PrefixJaccardScore.hpp>
#include <networkit/sparsification/SimmelianOverlapScore.hpp>
#include <networkit/sparsification/MultiscaleScore.hpp>
#include <networkit/sparsification/LocalSimilarityScore.hpp>
#include <networkit/sparsification/RandomEdgeScore.hpp>
#include <networkit/sparsification/GlobalThresholdFilter.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/auxiliary/Log.hpp>

namespace NetworKit {

class SparsificationBenchmark: public testing::Test {
protected:
    const int64_t n {250};

public:
    Graph makeCompleteGraph(count n) {
        Graph G(n);
        G.forNodePairs([&](node u, node v){
            G.addEdge(u,v);
        });
        G.shrinkToFit();
        return G;
    }
};

TEST_F(SparsificationBenchmark, completeGraphSimmelianSparsificationParametric) {
    int64_t n = this->n;
    Aux::Timer runtime;

    Graph G = SparsificationBenchmark::makeCompleteGraph(n);
    G.indexEdges();

    runtime.start();

    ChibaNishizekiTriangleEdgeScore counter(G);
    counter.run();
    std::vector<count> counts = counter.scores();

    SimmelianOverlapScore overlapScore(G, counts, 10);
    overlapScore.run();
    auto scores = overlapScore.scores();

    runtime.stop();
    INFO("[DONE] completeGraphSimmelianSparsificationParametric (" , runtime.elapsed().count() , " ms)");
}

TEST_F(SparsificationBenchmark, completeGraphSimmelianSparsificationNonParametric) {
    int64_t n = this->n;
    Aux::Timer runtime;

    Graph G = SparsificationBenchmark::makeCompleteGraph(n);
    G.indexEdges();

    runtime.start();

    ChibaNishizekiTriangleEdgeScore counter(G);
    counter.run();
    std::vector<count> counts = counter.scores();

    PrefixJaccardScore<count> jaccard(G, counts);
    jaccard.run();
    auto attribute = jaccard.scores();

    runtime.stop();
    INFO("[DONE] SimmelianSparsificationNonParametric (" , runtime.elapsed().count() , " ms)");
}

TEST_F(SparsificationBenchmark, completeGraphMultiscaleSparsification) {
    int64_t n = this->n;
    Aux::Timer runtime;

    Graph G = SparsificationBenchmark::makeCompleteGraph(n);
    G.indexEdges();

    runtime.start();

    std::vector<double> weight(G.upperEdgeIdBound());
    G.forEdges([&](node u, node v, edgeid eid) {
        weight[eid] = G.weight(u, v);
    });

    MultiscaleScore scorer(G, weight);
    scorer.run();
    auto scores = scorer.scores();

    runtime.stop();
    INFO("[DONE] MultiscaleSparsification (" , runtime.elapsed().count() , " ms)");
}

TEST_F(SparsificationBenchmark, completeGraphLocalSimilaritySparsification) {
    int64_t n = this->n;
    Aux::Timer runtime;

    Graph G = SparsificationBenchmark::makeCompleteGraph(n);
    G.indexEdges();

    runtime.start();

    ChibaNishizekiTriangleEdgeScore counter(G);
    counter.run();
    std::vector<count> triangles = counter.scores();

    LocalSimilarityScore localSimScore(G, triangles);
    localSimScore.run();
    auto attribute = localSimScore.scores();

    runtime.stop();
    INFO("[DONE] LocalSimilaritySparsification (" , runtime.elapsed().count() , " ms)");
}

TEST_F(SparsificationBenchmark, SparsificationBenchmarkGraphFile) {
    std::string path = "";

    std::cout << "[INPUT] .graph file path >" << std::endl;
    std::getline(std::cin, path);
    Aux::Timer runtime;

    // --------- IO
    std::cout << "[BEGIN] reading graph: " << path << std::endl;
    runtime.start();
    METISGraphReader reader;
    Graph g = reader.read(path);
    runtime.stop();
    std::cout << "[DONE] reading graph " << runtime.elapsedTag() << std::endl;

    // --------- Edge indexing
    std::cout << "[BEGIN] edge indexing: " << std::endl;
    runtime.start();
    g.indexEdges();
    runtime.stop();
    std::cout << "[DONE] edge indexing " << runtime.elapsedTag() << std::endl;

    // --------- Triangle counting
    std::cout << "[BEGIN] triangle counting: " << std::endl;
    runtime.start();
    ChibaNishizekiTriangleEdgeScore oldTriangleAttributizer(g);
    oldTriangleAttributizer.run();
    std::vector<count> oldTriangles = oldTriangleAttributizer.scores();
    runtime.stop();
    std::cout << "[DONE] Chiba Nishizeki triangle counting " << runtime.elapsedTag() << std::endl;

    // --------- Triangle counting
    std::cout << "[BEGIN] triangle counting: " << std::endl;
    runtime.start();
    TriangleEdgeScore triangleScore(g);
    triangleScore.run();
    std::vector<count> triangles = triangleScore.scores();
    runtime.stop();
    std::cout << "[DONE] triangle counting " << runtime.elapsedTag() << std::endl;

    // --------- Multiscale
    std::cout << "[BEGIN] multiscale attribute: " << std::endl;
    runtime.start();
    MultiscaleScore multiscaleScorer(g, std::vector<double>(triangles.begin(), triangles.end()));
    multiscaleScorer.run();
    std::vector<double> multiscale = multiscaleScorer.scores();
    runtime.stop();
    std::cout << "[DONE] multiscale attribute " << runtime.elapsedTag() << std::endl;

    std::cout << "[BEGIN] global filter (multiscale attribute): " << std::endl;
    runtime.start();
    GlobalThresholdFilter filter(g, multiscale, 0.5, false);
    Graph b = filter.calculate();
    runtime.stop();
    std::cout << "[DONE] global filter (multiscale attribute) " << runtime.elapsedTag() << std::endl;

    // --------- Simmelian Sparsification (Jaccard)
    std::cout << "[BEGIN] Simmelian Jaccard attribute: " << std::endl;
    runtime.start();
    PrefixJaccardScore<count> jaccardAttributizer(g, triangles);
    jaccardAttributizer.run();
    std::vector<double> jaccard = jaccardAttributizer.scores();
    runtime.stop();
    std::cout << "[DONE] Simmelian Jaccard attribute " << runtime.elapsedTag() << std::endl;

    std::cout << "[BEGIN] global filter (simmelian jaccard attribute): " << std::endl;
    runtime.start();
    GlobalThresholdFilter filter2(g, jaccard, 0.5, true);
    b = filter2.calculate();
    runtime.stop();
    std::cout << "[DONE] global filter (simmelian jaccard attribute) " << runtime.elapsedTag() << std::endl;

    // --------- Simmelian Sparsification (Overlap)
    std::cout << "[BEGIN] Simmelian Overlap attribute: " << std::endl;
    runtime.start();
    SimmelianOverlapScore overlapScore(g, triangles, 10);
    std::vector<double> overlap = overlapScore.scores();
    runtime.stop();
    std::cout << "[DONE] Simmelian Overlap attribute " << runtime.elapsedTag() << std::endl;

    std::cout << "[BEGIN] global filter (simmelian overlap attribute): " << std::endl;
    runtime.start();
    GlobalThresholdFilter filter3(g, overlap, 5, true);
    b = filter3.calculate();
    runtime.stop();
    std::cout << "[DONE] global filter (simmelian overlap attribute) " << runtime.elapsedTag() << std::endl;

    // --------- Local similarity Sparsification
    std::cout << "[BEGIN] Local Similarity attribute: " << std::endl;
    runtime.start();
    LocalSimilarityScore localSimScore(g, triangles);
    localSimScore.run();
    std::vector<double> minExponent = localSimScore.scores();
    runtime.stop();
    std::cout << "[DONE] Local Similarity attribute " << runtime.elapsedTag() << std::endl;

    std::cout << "[BEGIN] global filter (local similarity attribute): " << std::endl;
    runtime.start();
    GlobalThresholdFilter filter4 (g, minExponent, 0.37, true);
    b = filter4.calculate();
    runtime.stop();
    std::cout << "[DONE] global filter (local similarity attribute) " << runtime.elapsedTag() << std::endl;
}

} /* namespace NetworKit */
