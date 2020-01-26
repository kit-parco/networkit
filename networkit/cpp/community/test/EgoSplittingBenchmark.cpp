/*
 * EgoSplittingBenchmark.cpp
 *
 * Created: 2019-10-18
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/community/EgoSplitting.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class EgoSplittingBenchmark : public testing::Test {
public:
    void SetUp() {
        Aux::Random::setSeed(435913, false);
    }

    EgoSplittingBenchmark() {
        std::string graphPath;
        std::cout << "[INPUT] graph file path (edge list tab 0, like SNAP) > " << std::endl;
        std::getline(std::cin, graphPath);
        SNAPGraphReader reader;
        testGraph = reader.read(graphPath);
    }

    void benchEgoSplitting(const std::map<std::string, std::string> &parameters) {
        bool egoNetsParallel = true;
        EgoSplitting algo(testGraph, egoNetsParallel);
        algo.setParameters(parameters);

        algo.run();
        Cover cover = algo.getCover();

        std::cout << algo.timingsAsString() << std::endl;
    }

    Graph testGraph;
};

TEST_F(EgoSplittingBenchmark, benchNoExtend) {
    std::map<std::string, std::string> parameters;
    parameters["Extend EgoNet"] = "No";
    benchEgoSplitting(parameters);
}

TEST_F(EgoSplittingBenchmark, benchExtend) {
    std::map<std::string, std::string> parameters;
    parameters["Extend EgoNet"] = "Yes";
    benchEgoSplitting(parameters);
}


} /* namespace NetworKit */
