// no-networkit-format
/*
 * LAMGGTest.cpp
 *
 *  Created on: 20.11.2014
 *      Author: Michael
 */

#include <gtest/gtest.h>

#include <networkit/algebraic/Vector.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/generators/BarabasiAlbertGenerator.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/io/LineFileReader.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/METISGraphWriter.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/numerics/LAMG/MultiLevelSetup.hpp>
#include <networkit/numerics/LAMG/SolverLamg.hpp>
#include <networkit/numerics/GaussSeidelRelaxation.hpp>

namespace NetworKit {

class LAMGGTest : public testing::Test {
protected:
    const std::vector<std::string> GRAPH_INSTANCES = {"input/jazz.graph", "input/power.graph"};

    Vector randZeroSum(const Graph& graph, size_t seed) const;
    Vector randVector(count dimension) const;
};

TEST_F(LAMGGTest, testSmallGraphs) {
    METISGraphReader reader;
    GaussSeidelRelaxation<CSRMatrix> gaussSmoother, smoother;
    MultiLevelSetup<CSRMatrix> setup(gaussSmoother);
    Aux::Timer timer;
    for (index i = 0; i < GRAPH_INSTANCES.size(); ++i) {
        std::string graph = GRAPH_INSTANCES[i];
        Graph G = reader.read(graph);
        ConnectedComponents con(G);
        con.run();
        Partition comps = con.getPartition();
        if (comps.numberOfSubsets() > 1) { // disconnected graphs are currently not supported
            continue;
        }

        LevelHierarchy<CSRMatrix> hierarchy;
        timer.start();
        setup.setup(G, hierarchy);
        SolverLamg<CSRMatrix> solver(hierarchy, smoother);
        timer.stop();
        DEBUG("setup time\t ", timer.elapsedMilliseconds());

        Vector b(G.numberOfNodes());
        Vector x(G.numberOfNodes());

        b = randZeroSum(G, 12345);
        x = randVector(G.numberOfNodes());


        LAMGSolverStatus status;
        status.desiredResidualReduction = 1e-6 * b.length() / (hierarchy.at(0).getLaplacian() * x - b).length(); // needed for getting a relative residual <= 1e-6

        Vector result = x;
        DEBUG("Solving equation system - Gauss-Seidel");
        timer.start();
        solver.solve(result, b, status);
        timer.stop();

        EXPECT_TRUE(status.converged);

        DEBUG("solve time\t ", timer.elapsedMilliseconds());
        DEBUG("final residual = ", status.residual);
        DEBUG("numIters = ", status.numIters);
        DEBUG("DONE");

    }
}

Vector LAMGGTest::randVector(count dimension) const {
    Vector randVector(dimension);
    for (index i = 0; i < dimension; ++i) {
        randVector[i] = 2.0 * Aux::Random::probability() - 1.0;
    }

    // introduce bias
    for (index i = 0; i < dimension; ++i) {
        randVector[i] = randVector[i] * randVector[i];
    }

    return randVector;
}

Vector LAMGGTest::randZeroSum(const Graph& G, size_t seed) const {
    std::mt19937 rand(seed);
    auto rand_value = std::uniform_real_distribution<double>(-1.0, 1.0);
    ConnectedComponents con(G);
    count n = G.numberOfNodes();
    con.run();
    Partition comps = con.getPartition();

    /* Fill each component randomly such that its sum is 0 */
    Vector b(n, 0.0);

    for (int id : comps.getSubsetIds()) {
        auto indexes = comps.getMembers(id);
        assert(!indexes.empty());
        double sum = 0.0;
        for (auto entry : indexes) {
            b[entry] = rand_value(rand);
            sum += b[entry];
        }
        b[*indexes.begin()] -= sum;
    }

    return b;
}

} /* namespace NetworKit */
