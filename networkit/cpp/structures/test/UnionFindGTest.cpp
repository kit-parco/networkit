// no-networkit-format
/*
 * PartitionGTest.cpp
 *
 *  Created on: 04.12.2013
 *      Author: Maximilian Vogel (uocvf@student.kit.edu)
 */

#include <iostream>

#include <gtest/gtest.h>

#include <networkit/structures/UnionFind.hpp>


namespace NetworKit {

class UnionFindGTest: public testing::Test {};

TEST_F(UnionFindGTest, testAllToSingletons) {
    UnionFind p(10);
    p.allToSingletons();
    for (int i = 0; i < 10; i++) {
        EXPECT_TRUE(p.find(i) != p.find((i+1)%10));
    }
}

TEST_F(UnionFindGTest, testMergeSimple) {
    UnionFind p(10);
    p.merge(3,8);
    EXPECT_EQ(p.find(8),p.find(3));	
}


TEST_F(UnionFindGTest, testMergeSubsets) {
    UnionFind p(10);
    p.allToSingletons();
    p.merge(0,9);
    p.merge(1,8);
    p.merge(2,7);
    p.merge(0,1);
    p.merge(1,2);
    EXPECT_EQ(p.find(9),p.find(7));
}

TEST_F(UnionFindGTest, testMergeCircular) {
    UnionFind p(16);

    p.merge(0, 4);
    p.merge(1, 5);
    p.merge(2, 6);
    p.merge(3, 7);
    p.merge(8, 12);
    p.merge(9, 13);
    p.merge(10, 14);
    p.merge(11, 15);

    p.merge(0, 8);
    p.merge(1, 9);
    p.merge(2, 10);
    p.merge(3, 11);
    p.merge(4, 12);
    p.merge(5, 13);
    p.merge(6, 14);
    p.merge(7, 15);

    for (index i = 0; i < 15; ++i) {
        EXPECT_EQ(p.find(i), p.find((i+4) % 16));
        EXPECT_EQ(p.find(i), p.find((i+8) % 16));
        EXPECT_EQ(p.find(i), p.find((i+12) % 16));
        EXPECT_NE(p.find(i), p.find((i+1) % 16));
        EXPECT_NE(p.find(i), p.find((i+2) % 16));
        EXPECT_NE(p.find(i), p.find((i+3) % 16));
    }
}

} /* namespace NetworKit */
