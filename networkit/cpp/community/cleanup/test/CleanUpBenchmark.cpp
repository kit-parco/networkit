/*
 * CleanUpBenchmark.cpp
 *
 * Created: 2019-09-29
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/community/EgoSplitting.hpp>
#include <networkit/community/cleanup/SignificanceCommunityCleanUp.hpp>
#include <networkit/io/CoverReader.hpp>

namespace NetworKit {

class CleanUpBenchmark : public testing::Test {
public:
    void SetUp() {
        Aux::Random::setSeed(435913, false);
    }

    CleanUpBenchmark() {
        std::cout << "[INPUT] graph file path (edge list tab 0, like SNAP) > " << std::endl;
        std::getline(std::cin, graphsPath);
    }

    void benchCleanup(const std::string &graphName, bool mergeDiscarded) {
        const std::string filePrefix = graphsPath + "com-" + graphName;
        const std::string graphFile = filePrefix + ".ungraph.txt";
        const std::string coverFile = filePrefix + ".detected.cmty.txt";
        Aux::Log::setLogLevel("INFO");
        INFO("Read graph...");
        EdgeListReader reader('\t', 0);
        Graph G = reader.read(graphFile);
        INFO("Read cover...");
        CoverReader coverReader;
        Cover cover = coverReader.read(coverFile, G);
        std::vector<std::vector<node>> communities(cover.upperBound());
        cover.forEntries([&](node u, const std::set<index>& coms) {
            for (index s : coms) {
                communities[s].push_back(u);
            }
        });

        INFO("Start Cleanup...");

        Aux::Timer timer;
        timer.start();
        StochasticDistributionCalculator dist(2 * G.numberOfEdges() + G.numberOfNodes());
        SignificanceCommunityCleanUp cleanUp(G, communities, dist, 0.1, 0.1, 0.5, mergeDiscarded);
        cleanUp.run();
        timer.stop();
        std::cout << "Cleanup took " << timer.elapsedMilliseconds() << "ms" << std::endl;
    }

    std::string graphsPath;
};

TEST_F(CleanUpBenchmark, benchCommunityCleanup) {
    EdgeListReader reader(' ', 0);
    Graph G = reader.read("input/lfr_om3.graph");

    Aux::Timer timer;
    timer.start();
    EgoSplitting algo(G);
    algo.run();
    Cover cover = algo.getCover();
    std::vector<std::vector<node>> communities(cover.upperBound());
    cover.forEntries([&](node u, const std::set<index>& coms) {
        for (index s : coms) {
            communities[s].push_back(u);
        }
    });
    timer.stop();
    std::cout << "egosplitting took " << timer.elapsedMilliseconds() << "ms" << std::endl;

    timer.start();
    StochasticDistributionCalculator dist(2 * G.numberOfEdges() + G.numberOfNodes());
    SignificanceCommunityCleanUp cleanUp(G, communities, dist, 0.1, 0.1, 0.5);
    cleanUp.run();
    timer.stop();
    std::cout << "Cleanup took " << timer.elapsedMilliseconds() << "ms" << std::endl;

}

TEST_F(CleanUpBenchmark, benchLivejournal) {
    benchCleanup("lj", true);
}

TEST_F(CleanUpBenchmark, benchLivejournalNoMerge) {
    benchCleanup("lj", false);
}


} /* namespace NetworKit */

