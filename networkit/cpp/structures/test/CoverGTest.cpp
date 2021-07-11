// no-networkit-format
/*
 * CoverGTest.cpp
 *
 *  Created on: 12.12.2013
 *      Author: Maximilian Vogel (uocvf@student.kit.edu)
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/structures/Cover.hpp>

#include <iostream>

namespace NetworKit {

class CoverGTest: public testing::Test {};

TEST_F(CoverGTest, testConstructor) {
    Cover c(10);
    EXPECT_EQ(0u, c.lowerBound());
    EXPECT_EQ(1u, c.upperBound());
}

TEST_F(CoverGTest, testAllToSingletonsAndUpperBound) {
    Cover c(10);
    EXPECT_EQ(1u, c.upperBound());
    c.allToSingletons();
    EXPECT_EQ(0u, c.lowerBound());
    EXPECT_EQ(11u, c.upperBound());
}

TEST_F(CoverGTest, testContains) {
    Cover c(10);
    c.toSingleton(0);
    EXPECT_TRUE(c.contains(0));
    EXPECT_FALSE(c.contains(1));
}

TEST_F(CoverGTest, testUpperBoundAfterMerges) {
    count n = 10;
    Cover c(n);
    for (index i = 0; i < n; i+=2) {
        index sid = c.toSingleton(i);
        c.toSingleton(i+1);
        c.addToSubset(sid,i+1);
    }
    for (index i = 0; i < n; i++) {
        c.addToSubset(i+1,0);
    }
    c.mergeSubsets(1,3);
    c.mergeSubsets(5,11);
    EXPECT_EQ(13u, c.upperBound());
}


TEST_F(CoverGTest, testToSingleton) {
    count n = 10;
    Cover c(n);
    c.allToSingletons();
    std::set<index> controlSet2;
    controlSet2.insert(1);
    DEBUG("c[0] ", c[0], " and controlSet2 ", controlSet2);
    EXPECT_TRUE(c[0] == controlSet2);
    c.addToSubset(5,0);
    c.addToSubset(2,0);
    c.addToSubset(3,0);
    c.addToSubset(4,0);
    c.addToSubset(0,1);
    c.toSingleton(0);
    std::set<index> controlSet;
    controlSet.insert(11);
    DEBUG("c[0] ", c[0], " and controlSet ", controlSet);
    EXPECT_TRUE(c[0] == controlSet);
}

TEST_F(CoverGTest, testAddToSubset) {
    count n = 10;
    Cover c(n);
    c.addToSubset(0,0);
    c.addToSubset(0,1);
    std::set<index> controlSet = {0};
    EXPECT_TRUE(c.inSameSubset(0,1));
    EXPECT_TRUE(c[0] == controlSet);
    EXPECT_TRUE(c[1] == controlSet);
}


TEST_F(CoverGTest, testAddToSubset2) {
    count n = 10;
    Cover c(n);
    index sid = c.toSingleton(0);
    c.addToSubset(sid,5);
    EXPECT_TRUE(c.inSameSubset(0,5));
}

TEST_F(CoverGTest, testMoveToSubset) {
    count n = 10;
    Cover c(n);
    c.allToSingletons();
    c.addToSubset(5,0);
    c.addToSubset(2,0);
    c.addToSubset(3,0);
    c.addToSubset(4,0);
    c.addToSubset(0,1);
    c.moveToSubset(8,0);
    std::set<index> controlSet = {8};
    EXPECT_EQ(c[0],controlSet);
}

TEST_F(CoverGTest, testSubsetSizesWithUnassignedElements) {
    count n = 10;
    Cover c(n);
    for (index i = 0; i < n; i+=2) {
        c.toSingleton(i);
    }
    std::vector<index> controlSet = {1,1,1,1,1};
    EXPECT_EQ(c.subsetSizes(),controlSet);
}

TEST_F(CoverGTest, testSubsetSizesTrivial) {
    count n = 10;
    Cover c(n);
    for (index i = 0; i < n; i++) {
        c.toSingleton(i);
    }
    std::vector<index> controlSet = {1,1,1,1,1,1,1,1,1,1};
    EXPECT_EQ(c.subsetSizes(),controlSet);
}

TEST_F(CoverGTest, testSubsetSizesTrivial2) {
    count n = 10;
    Cover c(n);
    for (index i = 0; i < n; i+=2) {
        c.toSingleton(i);
    }
    for (index i = 1; i < n; i+=2) {
        c.addToSubset((i/2)+1,i);
    }
    std::vector<index> controlSet = {2,2,2,2,2};
    EXPECT_EQ(c.subsetSizes(),controlSet);
}

TEST_F(CoverGTest, testSubsetSizesAssignedToMultipleSubsets) {
    count n = 10;
    Cover c(n);
    for (index i = 0; i < n; i++) {
        c.toSingleton(i);
    }
    for (index i = 1; i < n; i+=2) {
        c.addToSubset(i,i);
    }
    std::vector<index> controlSet = {2,1,2,1,2,1,2,1,2,1};
    EXPECT_EQ(c.subsetSizes(),controlSet);
}

TEST_F(CoverGTest, testSubsetSizesAssignedToMultipleSubsets2) {
    count n = 10;
    Cover c(n);
    for (index i = 0; i < n; i+=2) {
        index sid = c.toSingleton(i);
        c.toSingleton(i+1);
        c.addToSubset(sid,i+1);
    }
    for (index i = 0; i < n; i++) {
        c.addToSubset(i+1,0);
    }
    std::vector<index> controlSet = {2,2,3,2,3,2,3,2,3,2};
    EXPECT_EQ(c.subsetSizes(),controlSet);
}

TEST_F(CoverGTest, testSubsetSizeMapMultipleSets) {
    count n = 10;
    Cover c(n);
    for (index i = 0; i < n; i++) {
        c.toSingleton(i);
    }
    for (index i = 1; i < n; i+=2) {
        c.addToSubset(i,i);
    }
    std::map<index,count> controlMap;
    controlMap[1] = 2;
    controlMap[2] = 1;
    controlMap[3] = 2;
    controlMap[4] = 1;
    controlMap[5] = 2;
    controlMap[6] = 1;
    controlMap[7] = 2;
    controlMap[8] = 1;
    controlMap[9] = 2;
    controlMap[10] = 1;
    EXPECT_EQ(c.subsetSizeMap(),controlMap);
}

TEST_F(CoverGTest, testMergeSubsetsAndGetMembers) {
    count n = 10;
    Cover c(n);
    for (index i = 0; i < n; i+=2) {
        index sid = c.toSingleton(i);
        c.toSingleton(i+1);
        c.addToSubset(sid,i+1);
    }
    for (index i = 0; i < n; i++) {
        c.addToSubset(i+1,0);
    }
    c.mergeSubsets(1,3);
    c.mergeSubsets(5,11);
    auto c11 = c.getMembers(11);
    std::vector<index> controlSetSizes = {2,2,2,3,2,3,2,6};
    // remaining subset IDs 2,4,6,7,8,9,10,12
    // their sizes          2,2,2,3,2,3,2,6
    std::set<index> controlSetMembers = {0,1,2,3,4,5};
    EXPECT_EQ(controlSetSizes,c.subsetSizes()); // check if subsets sizes are correct
    auto c12 = c.getMembers(12);
    EXPECT_EQ(c12,controlSetMembers); // check if elements of merged subsets are correct
}

TEST_F(CoverGTest, testNumberOfSubsets) {
    count n = 10;
    Cover c(n);
    for (index i = 0; i < n; i+=2) {
        index sid = c.toSingleton(i);
        c.toSingleton(i+1);
        c.addToSubset(sid,i+1);
    }
    for (index i = 0; i < n; i++) {
        c.addToSubset(i+1,0);
    }
    EXPECT_EQ(n,c.numberOfSubsets());
    c.mergeSubsets(1,2);
    c.mergeSubsets(3,11);
    EXPECT_EQ(8u, c.numberOfSubsets());
}


TEST_F(CoverGTest, testSubsetsOf) {
    count n = 10;
    Cover c(n);
    for (index i = 0; i < n; i+=2) {
        index sid = c.toSingleton(i);
        c.toSingleton(i+1);
        c.addToSubset(sid,i+1);
    }
    for (index i = 0; i < n; i++) {
        c.addToSubset(i+1,0);
    }
    auto subsetsOf0 = c.subsetsOf(0);
    std::set<index> controlSet0 = {1,2,3,4,5,6,7,8,9,10};
    auto subsetsOf3 = c.subsetsOf(3);
    std::set<index> controlSet3 = {3,4};
    EXPECT_EQ(controlSet0,subsetsOf0);
    EXPECT_EQ(controlSet3,subsetsOf3);
    c.mergeSubsets(1,3);
    c.mergeSubsets(5,11);
    subsetsOf0 = c.subsetsOf(0);
    controlSet0 = {2,4,6,7,8,9,10,12};
    EXPECT_EQ(controlSet0,subsetsOf0);
}

TEST_F(CoverGTest, testInSameSubset) {
    count n = 10;
    Cover c(n);
    for (index i = 0; i < n; i+=2) {
        index sid = c.toSingleton(i);
        c.toSingleton(i+1);
        c.addToSubset(sid,i+1);
    }
    EXPECT_TRUE(c.inSameSubset(0,1));
    EXPECT_FALSE(c.inSameSubset(0,2));
    EXPECT_FALSE(c.inSameSubset(1,5));
    c.mergeSubsets(1,3);
    EXPECT_TRUE(c.inSameSubset(0,1));
    EXPECT_TRUE(c.inSameSubset(0,2));
    EXPECT_FALSE(c.inSameSubset(1,5));
    c.mergeSubsets(5,11);
    EXPECT_TRUE(c.inSameSubset(0,1));
    EXPECT_TRUE(c.inSameSubset(0,2));
    EXPECT_TRUE(c.inSameSubset(1,5));
}


} /* namespace NetworKit */
