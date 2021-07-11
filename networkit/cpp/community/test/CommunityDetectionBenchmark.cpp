// no-networkit-format
/*
 * CommunityDetectionBenchmark.h
 *
 *  Created on: 16.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <gtest/gtest.h>

#include <map>
#include <functional>

#include <networkit/community/PLP.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/centrality/Betweenness.hpp>
#include <networkit/centrality/PageRank.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/structures/Partition.hpp>

#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class CommunityDetectionBenchmark: public testing::Test {
public:
    virtual ~CommunityDetectionBenchmark() = default;

protected:
    METISGraphReader metisReader;

};

constexpr int runs = 12;

TEST_F(CommunityDetectionBenchmark, benchClusteringAlgos) {
    Aux::Timer timer;
    Modularity mod;

    // std::string graph = "../graphs/in-2004.graph";
    // std::string graph = "../graphs/uk-2002.graph";
    // std::string graph = "../graphs/uk-2007-05.graph";
    std::string graph = "input/polblogs.graph";

    DEBUG("Reading graph file ", graph.c_str(), " ...");
    timer.start();
    const Graph G = this->metisReader.read(graph);
    timer.stop();
    DEBUG("Reading graph took ", timer.elapsedMilliseconds() / 1000.0, "s");

    for (int r = 0; r < runs; r++) {
        Graph Gcopy = G;
        PLP algo(Gcopy);

        timer.start();
        algo.run();
        Partition zeta = algo.getPartition();
        timer.stop();

        auto communitySizes = zeta.subsetSizes();


        INFO("Parallel Label Propagation on ", graph.c_str(), ": ",
                (timer.elapsedMilliseconds() / 1000.0),
                "s,\t#communities: ", zeta.numberOfSubsets(),
                ",\tmodularity: ", mod.getQuality(zeta, G));
    }

    for (int r = 0; r < runs; r++) {
        Graph Gcopy = G;
        PLM algo(Gcopy);

        timer.start();
        algo.run();
        Partition zeta = algo.getPartition();
        timer.stop();

        auto communitySizes = zeta.subsetSizes();

        INFO("Parallel Louvain on ", graph.c_str(), ": ",
                (timer.elapsedMilliseconds() / 1000.0),
                "s,\t#communities: ", zeta.numberOfSubsets(),
                ",\tmodularity: ", mod.getQuality(zeta, G));
    }
}

TEST_F(CommunityDetectionBenchmark, benchPageRankCentrality) {
    Aux::Timer timer;

    // std::string graph = "../graphs/uk-2002.graph";
    std::string graph = "input/polblogs.graph";

const Graph G = this->metisReader.read(graph);

    for (int r = 0; r < runs; r++) {
        PageRank cen(G, 1e-6);

        timer.start();
        cen.run();
        timer.stop();
        auto ranking = cen.ranking();

        INFO("Page Rank Centrality on ", graph.c_str(), ": ",
                (timer.elapsedMilliseconds() / 1000.0),
                "s,\t ranking: [(",
                ranking[0].first, ": ", ranking[0].second, "), (",
                ranking[1].first, ": ", ranking[1].second, ") ...]");
    }
}

TEST_F(CommunityDetectionBenchmark, benchBetweennessCentrality) {
    Aux::Timer timer;

    // std::string graph = "../graphs/cond-mat-2005.graph";
    std::string graph = "input/polblogs.graph";

const Graph G = this->metisReader.read(graph);

    for (int r = 0; r < runs; r++) {
        Betweenness cen(G);

        timer.start();
        cen.run();
        timer.stop();
        auto ranking = cen.ranking();

        INFO("Betweenness Centrality on ", graph.c_str(), ": ",
                (timer.elapsedMilliseconds() / 1000.0),
                "s,\t ranking: [(",
                ranking[0].first, ": ", ranking[0].second, "), (",
                ranking[1].first, ": ", ranking[1].second, ") ...]");

    }
}


} /* namespace NetworKit */

