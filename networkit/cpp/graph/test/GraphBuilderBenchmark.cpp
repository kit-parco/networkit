// no-networkit-format
/*
 * GraphBuilderBenchmark.cpp
 *
 *  Created on: 04.12.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com), Manuel Penschuck <networkit@manuel.jetzt>
 *
 */

#ifndef NOGTEST

#include <networkit/graph/test/GraphBuilderBenchmark.hpp>
#include <networkit/graph/GraphBuilder.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/generators/ErdosRenyiEnumerator.hpp>

namespace NetworKit {

constexpr node graphBuilderNodes = 200000;
constexpr double graphBuilderProb = 0.0001;

TEST_F(GraphBuilderBenchmark, benchmarkGraphBuilderBaseline) {
    ErdosRenyiEnumerator<> ere(graphBuilderNodes, graphBuilderProb, false);
    auto t1 = timeOnce([&]() {
        ere.forEdgesParallel([&](int tid, node u, node v) {volatile auto tmp = u;});
    });

    std::cout << "ErdosRenyi only:\t\t" << t1 << " ms\n";
}

TEST_F(GraphBuilderBenchmark, benchmarkGraphBuilderParFillSeqBuild) {
    GraphBuilder builder(graphBuilderNodes);
    count m_actual = 0;
    ErdosRenyiEnumerator<> ere(graphBuilderNodes, graphBuilderProb, false);
    auto t1 = timeOnce([&]() {
        ere.forEdgesParallel([&](int tid, node u, node v) {
            builder.addHalfEdge(u, v);
        });
    });
    auto t2 = timeOnce([&]() {
        auto G = builder.toGraph(true, false);
        m_actual = G.numberOfEdges();
    });
    EXPECT_NEAR(m_actual, ere.expectedNumberOfEdges(), 0.1 * ere.expectedNumberOfEdges());
    std::cout << "parallelForNodePairs + toGraphSequentiel:\t\t" << t1 << " + " << t2 << " = " << (t1 + t2) << " ms\n";
}

TEST_F(GraphBuilderBenchmark, benchmarkGraphBuilderParFillParBuild) {
    GraphBuilder builder(graphBuilderNodes);
    count m_actual = 0;
    ErdosRenyiEnumerator<> ere(graphBuilderNodes, graphBuilderProb, false);
    auto t1 = timeOnce([&]() {
        ere.forEdgesParallel([&](int tid, node u, node v) {
            builder.addHalfEdge(u, v);
        });
    });
    auto t2 = timeOnce([&]() {
        auto G = builder.toGraph(true, true);
        m_actual = G.numberOfEdges();
    });
    EXPECT_NEAR(m_actual, ere.expectedNumberOfEdges(), 0.1 * ere.expectedNumberOfEdges());
    std::cout << "parallelForNodePairs + toGraphParallel:\t\t" << t1 << " + " << t2 << " = " << (t1 + t2) << " ms\n";
}


TEST_F(GraphBuilderBenchmark, benchmarkMETISReader) {
    METISGraphReader reader;
    measureInMs([&]() {
        auto G = reader.read("../algoDaten/graphs/eu-2005.graph");
        return G.numberOfNodes();
    }, 20);
}

} /* namespace NetworKit */

#endif /*NOGTEST */
