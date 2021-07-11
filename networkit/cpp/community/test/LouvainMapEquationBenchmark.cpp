/*
 * LouvainMapEquationBenchmark.cpp
 *
 * Created on: 2019-10-31
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/community/LouvainMapEquation.hpp>
#include <networkit/generators/ClusteredRandomGraphGenerator.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>

namespace NetworKit {

class MapEquationBenchmark : public testing::Test {
public:
    void SetUp() { Aux::Random::setSeed(435913, false); }
};

TEST_F(MapEquationBenchmark, benchLarge) {
    ClusteredRandomGraphGenerator generator(5000, 200, 0.5, 0.002);
    Graph G = generator.generate();
    Partition groundTruth = generator.getCommunities();
    Aux::Timer timer{};
    timer.start();

    LouvainMapEquation mapequation(G, false);
    mapequation.run();
    auto partition = mapequation.getPartition();

    timer.stop();
    Aux::Log::setLogLevel("INFO");
    INFO(mapequation.toString(), " took ", timer.elapsedMilliseconds(), "ms");
}

TEST_F(MapEquationBenchmark, benchLargeHierarchical) {
    std::cout << "[INPUT] graph file path (edge list tab 0, like SNAP) > " << std::endl;
    std::string graphPath;
    std::getline(std::cin, graphPath);
    SNAPGraphReader reader;
    Graph G = reader.read(graphPath);
    Aux::Random::setSeed(420, true);
    Aux::Timer timer{};
    timer.start();

    LouvainMapEquation mapequation(G, true);
    mapequation.run();
    auto partition = mapequation.getPartition();

    timer.stop();
    INFO(mapequation.toString(), " took ", timer.elapsedMilliseconds(), "ms");
}

} // namespace NetworKit
