// no-networkit-format
/*
 * BiconnectedComponentsGTest.cpp
 *
 * Created on: March 2018
 *     Author: Eugenio Angriman
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/components/BiconnectedComponents.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class BiconnectedComponentsGTest : public testing::Test {};

TEST_F(BiconnectedComponentsGTest, testBiconnectedComponentsTiny) {
    Graph G(9, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(1, 4);
    G.addEdge(0, 5);
    G.addEdge(0, 6);
    G.addEdge(4, 5);
    G.addEdge(2, 3);
    G.addEdge(6, 8);
    G.addEdge(6, 7);
    G.addEdge(7, 8);
    BiconnectedComponents bc(G);
    bc.run();

    EXPECT_EQ(bc.numberOfComponents(), 4);
}

TEST_F(BiconnectedComponentsGTest, testBiconnectedComponents) {
    Aux::Random::setSeed(42, false);
    Graph G = ErdosRenyiGenerator(200, 0.01, false).generate();

    BiconnectedComponents bc(G);
    bc.run();

    for (auto component : bc.getComponents()) {
        std::unordered_set<node> subgraph(component.begin(), component.end());
        const auto G1 = GraphTools::subgraphFromNodes(G, subgraph);

        G1.forNodes([&](node v) {
            auto G2(G1);
            G2.removeNode(v);
            ConnectedComponents cc(G2);
            cc.run();
            EXPECT_EQ(cc.numberOfComponents(), 1);
        });
    }
}

} // namespace NetworKit
