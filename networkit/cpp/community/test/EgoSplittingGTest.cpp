/*
 * EgoSplittingGTest.cpp
 *
 * Created: 2019-10-15
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>
#include <functional>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/community/EgoSplitting.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class EgoSplittingGTest : public testing::Test {
public:
    void SetUp() {
        Aux::Random::setSeed(435913, false);
    }

    EgoSplittingGTest() {
        std::string inputDir = "input";
        METISGraphReader reader{};
        testGraph = reader.read(inputDir + "/lfr_small.graph");
        testGraph.removeNode(0);
        testGraph.forNeighborsOf(1, [&](node v) {
            testGraph.removeEdge(1, v);
        });
        testGraph.addEdge(1, 2);
    }

    void testEgoSplitting(const std::map<std::string, std::string> &parameters, bool parallelEgoNets = false) {
        EgoSplitting algo(testGraph, parallelEgoNets);
        algo.setParameters(parameters);
        algo.run();
        Cover cover = algo.getCover();

        EXPECT_GE(cover.numberOfSubsets(), 5);
        for (auto size : cover.subsetSizes()) {
            EXPECT_GT(size, 4) << "discard communities with 4 or less nodes";
            EXPECT_LT(size, 40);
        }
    }

    Graph testGraph;
};

TEST_F(EgoSplittingGTest, testEgoSplittingSimple) {
    std::map<std::string, std::string> parameters;
    parameters["Extend EgoNet"] = "No";
    parameters["Connect Personas Internally"] = "No";
    parameters["Cleanup"] = "No";
    testEgoSplitting(parameters);
}

TEST_F(EgoSplittingGTest, testEgoSplittingExtend) {
    std::map<std::string, std::string> parameters;
    parameters["Extend EgoNet"] = "Yes";
    parameters["Connect Personas Internally"] = "No";
    parameters["Cleanup"] = "No";
    testEgoSplitting(parameters);
}

TEST_F(EgoSplittingGTest, testEgoSplittingConnectPersonas) {
    std::map<std::string, std::string> parameters;
    parameters["Extend EgoNet"] = "Yes";
    parameters["Connect Personas Internally"] = "Yes";
    parameters["Cleanup"] = "No";
    testEgoSplitting(parameters);
}

TEST_F(EgoSplittingGTest, testEgoSplittingFull) {
    std::map<std::string, std::string> parameters;
    parameters["Extend EgoNet"] = "Yes";
    parameters["Connect Personas Internally"] = "Yes";
    parameters["Cleanup"] = "Yes";
    testEgoSplitting(parameters);
}

TEST_F(EgoSplittingGTest, testEgoSplittingParallel) {
    std::map<std::string, std::string> parameters;
    parameters["Extend EgoNet"] = "Yes";
    parameters["Connect Personas Internally"] = "Yes";
    parameters["Cleanup"] = "Yes";
    testEgoSplitting(parameters, true);
}

} /* namespace NetworKit */
