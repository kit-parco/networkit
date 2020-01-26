/*
 * NodeMappingGTest.cpp
 *
 * Created: 2020-01-12
 * Author: Armin Wiebigke
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include <networkit/structures/NodeMapping.hpp>

namespace NetworKit {

class NodeMappingGTest : public testing::Test {
};

TEST_F(NodeMappingGTest, testSimple) {
    NodeMapping nodeMapping(10);

    nodeMapping.addNode(2);
    nodeMapping.addNode(9);
    nodeMapping.addNode(5);

    ASSERT_EQ(nodeMapping.nodeCount(), 3);
    ASSERT_TRUE(nodeMapping.isMapped(2));
    ASSERT_TRUE(nodeMapping.isMapped(5));
    ASSERT_TRUE(nodeMapping.isMapped(9));

    std::vector<node> sortedNodes = nodeMapping.globalNodes();
    std::sort(sortedNodes.begin(), sortedNodes.end());
    ASSERT_EQ(sortedNodes, std::vector<node>({2, 5, 9}));

    ASSERT_EQ(nodeMapping.toLocal(2), 0);
    ASSERT_EQ(nodeMapping.toLocal(9), 1);
    ASSERT_EQ(nodeMapping.toLocal(5), 2);

    ASSERT_EQ(nodeMapping.toGlobal(0), 2);
    ASSERT_EQ(nodeMapping.toGlobal(1), 9);
    ASSERT_EQ(nodeMapping.toGlobal(2), 5);
}

TEST_F(NodeMappingGTest, testReset) {
    NodeMapping nodeMapping(10);

    nodeMapping.addNode(2);
    nodeMapping.addNode(9);
    nodeMapping.addNode(5);
    nodeMapping.reset();

    EXPECT_EQ(nodeMapping.nodeCount(), 0);
    EXPECT_FALSE(nodeMapping.isMapped(2));
    EXPECT_FALSE(nodeMapping.isMapped(5));
    EXPECT_FALSE(nodeMapping.isMapped(9));
}

TEST_F(NodeMappingGTest, testResetToIndex) {
    NodeMapping nodeMapping(10);

    nodeMapping.addNode(2);
    nodeMapping.addNode(9);
    nodeMapping.addNode(5);
    nodeMapping.addNode(3);
    nodeMapping.reset(2);

    EXPECT_EQ(nodeMapping.nodeCount(), 2);

    EXPECT_TRUE(nodeMapping.isMapped(2));
    EXPECT_TRUE(nodeMapping.isMapped(9));
    EXPECT_FALSE(nodeMapping.isMapped(5));
    EXPECT_FALSE(nodeMapping.isMapped(3));

    EXPECT_EQ(nodeMapping.toLocal(2), 0);
    EXPECT_EQ(nodeMapping.toLocal(9), 1);
    EXPECT_EQ(nodeMapping.toGlobal(0), 2);
    EXPECT_EQ(nodeMapping.toGlobal(1), 9);
}

}