// no-networkit-format
/*
 * PostscriptWriterGTest.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#include <gtest/gtest.h>

#include <vector>

#include <networkit/viz/PostscriptWriter.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/generators/ClusteredRandomGraphGenerator.hpp>
#include <networkit/io/PartitionWriter.hpp>
#include <networkit/io/DibapGraphReader.hpp>
#include <networkit/generators/PubWebGenerator.hpp>
#include <networkit/auxiliary/Random.hpp>

namespace NetworKit {

class VizGTest : public testing::Test {};


TEST_F(VizGTest, testPostscriptWriterOnRandomGraph) {
    // create graph
    count n = 60;
    count numClusters = 3;
    double pin = 0.35;
    double pout = 0.05;

    ClusteredRandomGraphGenerator graphGen(n, numClusters, pin, pout);
    Graph G = graphGen.generate();
    std::vector<Point2D> coordinates(G.upperNodeIdBound());

    // create coordinates
    G.forNodes([&](node u) {
        coordinates[u] = {Aux::Random::probability(), Aux::Random::probability()};
    });


    // write graph to file
    std::string path = "output/testGraph.eps";
    PostscriptWriter psWriter;
    psWriter.write(G, coordinates, path);

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

#ifndef NETWORKIT_WINDOWS
TEST_F(VizGTest, testPostscriptWriterOnRealGraph) {
    // read graph and coordinates from binary file
    DibapGraphReader reader;
    Graph G = reader.read("input/airfoil1.gi");
    const auto coordinates = reader.moveCoordinates();


    // write graph to file
    std::string path = "output/airfoil1.eps";
    PostscriptWriter psWriter;
    psWriter.write(G, Point<>::pointVectorToPoint2D(coordinates), path);

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}
#endif

} /* namespace NetworKit */
