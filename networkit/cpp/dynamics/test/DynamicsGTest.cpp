// no-networkit-format
/*
 * DynamicsGTest.cpp
 *
 *  Created on: 24.12.2013
 *      Author: cls
 */

#include <gtest/gtest.h>

#include <networkit/dynamics/DGSStreamParser.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/dynamics/GraphUpdater.hpp>
#include <networkit/dynamics/GraphDifference.hpp>

namespace NetworKit {

class DynamicsGTest: public testing::Test {};

TEST_F(DynamicsGTest, testDGSStreamParser) {
    DGSStreamParser parser("input/example2.dgs");
    auto stream = parser.getStream();
    ASSERT_EQ(stream.size(),16);
    ASSERT_EQ(stream[0].type,GraphEvent::NODE_ADDITION);
    ASSERT_EQ(stream[1].type,GraphEvent::NODE_ADDITION);
    ASSERT_EQ(stream[2].type,GraphEvent::EDGE_ADDITION);
    ASSERT_EQ(stream[3].type,GraphEvent::TIME_STEP);
    ASSERT_EQ(stream[4].type,GraphEvent::EDGE_WEIGHT_UPDATE);
    ASSERT_EQ(stream[5].type,GraphEvent::EDGE_REMOVAL);
    ASSERT_EQ(stream[6].type,GraphEvent::NODE_REMOVAL);
    ASSERT_EQ(stream[7].type,GraphEvent::NODE_REMOVAL);
    ASSERT_EQ(stream[8].type,GraphEvent::NODE_ADDITION);
    ASSERT_EQ(stream[9].type,GraphEvent::NODE_ADDITION);
    ASSERT_EQ(stream[10].type,GraphEvent::EDGE_ADDITION);
    ASSERT_EQ(stream[11].type,GraphEvent::NODE_ADDITION);
    ASSERT_EQ(stream[12].type,GraphEvent::EDGE_ADDITION);
    ASSERT_EQ(stream[13].type,GraphEvent::NODE_ADDITION);
    ASSERT_EQ(stream[14].type,GraphEvent::NODE_REMOVAL);
    ASSERT_EQ(stream[15].type,GraphEvent::NODE_RESTORATION);

    //apply updates
    Graph G(0,true);
    GraphUpdater updater(G);
    updater.update(stream);
    ASSERT_EQ(G.numberOfNodes(),4);
    ASSERT_EQ(G.numberOfEdges(),2);
}


TEST_F(DynamicsGTest, debugDGSStreamParserOnRealGraph) {
    std::string path;
    std::cout << "enter .dgs file path: ";
    std::cin >> path;
    DGSStreamParser parser(path);
    auto stream = parser.getStream();
}

TEST_F(DynamicsGTest, testGraphEventIncrement) {
    Graph G(2, true, false); //undirected
    Graph H(2, true, true); //directed
    G.addEdge(0, 1, 3.14);
    H.addEdge(0, 1, 3.14);
    GraphEvent event(GraphEvent::EDGE_WEIGHT_INCREMENT, 0, 1, 2.1);
    std::vector<GraphEvent> eventstream(1);
    eventstream.push_back(event);
    GraphUpdater Gupdater(G);
    GraphUpdater Hupdater(H);
    Gupdater.update(eventstream);
    Hupdater.update(eventstream);
    EXPECT_EQ(G.weight(0,1), 5.24);
    EXPECT_EQ(H.weight(0,1), 5.24);



}

namespace {
    // helper methods
    std::string edits_to_string(const std::vector<GraphEvent>& events) {
        std::stringstream ss;
        ss << "[";
        for (size_t i = 0; i < events.size(); ++i) {
            ss << events[i].toString();
            if (i < events.size() - 1) {
                ss << ", ";
            }
        }
        ss << "]";
        return ss.str();
    }

    void expect_graph_equals(const Graph& G1, const Graph &G2) {
        EXPECT_EQ(G1.numberOfNodes(), G2.numberOfNodes());
        EXPECT_EQ(G1.numberOfEdges(), G2.numberOfEdges());

        GraphDifference diff1(G1, G2);
        diff1.run();
        EXPECT_EQ(diff1.getNumberOfEdits(), 0);
        EXPECT_TRUE(diff1.getEdits().empty()) << edits_to_string(diff1.getEdits());

        GraphDifference diff2(G2, G2);
        diff2.run();
        EXPECT_EQ(diff2.getNumberOfEdits(), 0);
        EXPECT_TRUE(diff2.getEdits().empty()) << edits_to_string(diff2.getEdits());
    }
}

TEST_F(DynamicsGTest, testGraphDifference) {
    Graph G1(11, false, false);
    Graph G2(8, false, false);

    G1.addEdge(2, 4);
    G1.addEdge(2, 6);
    G1.addEdge(9, 10);
    G1.addEdge(4, 10);
    G1.removeNode(3);
    G1.removeNode(8);

    G2.addEdge(3, 2);
    G2.removeNode(4);

    {
        GraphDifference diff(G1, G2);
        diff.run();

        EXPECT_EQ(diff.getNumberOfEdgeRemovals(), 4) << edits_to_string(diff.getEdits());
        EXPECT_EQ(diff.getNumberOfNodeRemovals(), 3) << edits_to_string(diff.getEdits());
        EXPECT_EQ(diff.getNumberOfNodeRestorations(), 1) << edits_to_string(diff.getEdits());
        EXPECT_EQ(diff.getNumberOfEdgeAdditions(), 1) << edits_to_string(diff.getEdits());

        Graph H = G1;
        GraphUpdater Hupdater(H);
        Hupdater.update(diff.getEdits());
        EXPECT_TRUE(H.hasEdge(2, 3));

        expect_graph_equals(H, G2);
    }

    {
        GraphDifference diff(G2, G1);
        diff.run();

        EXPECT_EQ(diff.getNumberOfEdgeAdditions(), 4) << edits_to_string(diff.getEdits());
        EXPECT_EQ(diff.getNumberOfNodeRemovals(), 1) << edits_to_string(diff.getEdits());
        EXPECT_EQ(diff.getNumberOfNodeAdditions(), 2) << edits_to_string(diff.getEdits());
        EXPECT_EQ(diff.getNumberOfNodeRestorations(), 1) << edits_to_string(diff.getEdits());
        EXPECT_EQ(diff.getNumberOfEdgeRemovals(), 1) << edits_to_string(diff.getEdits());

        Graph H = G2;
        GraphUpdater Hupdater(H);
        Hupdater.update(diff.getEdits());

        expect_graph_equals(H, G1);
    }
}

TEST_F(DynamicsGTest, testGraphDifferenceDirectedWeighted) {
    Graph G1(8, true, true);
    Graph G2(8, true, true);

    G1.addEdge(0, 1, 0.5);
    G1.addEdge(1, 0, 2);
    G1.addEdge(2, 3, 2.0);

    G2.addEdge(1, 0, 0.5);
    G2.addEdge(3, 2, 1.0);
    G2.addEdge(3, 4, 2.5);

    for (const auto& gs : std::vector<std::pair<Graph, Graph>>({{G1, G2}, {G2, G1}})) {
        const Graph& g1 = gs.first;
        const Graph& g2 = gs.second;

        GraphDifference diff(g1, g2);
        diff.run();

        EXPECT_EQ(diff.getNumberOfEdgeRemovals(), 2) << edits_to_string(diff.getEdits());
        EXPECT_EQ(diff.getNumberOfEdgeAdditions(), 2) << edits_to_string(diff.getEdits());
        EXPECT_EQ(diff.getNumberOfEdgeWeightUpdates(), 1) << edits_to_string(diff.getEdits());
        EXPECT_EQ(diff.getNumberOfEdits(), 5) << edits_to_string(diff.getEdits());

        Graph H = g1;
        GraphUpdater Hupdater(H);
        Hupdater.update(diff.getEdits());

        expect_graph_equals(H, g2);
    }
}

} /* namespace NetworKit */
