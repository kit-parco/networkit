// no-networkit-format
/*
 * AuxGTest.cpp
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

// this define is an obscure fix for std::this_thread::sleep_for to work - the issue is described here: http://stackoverflow.com/questions/4438084/stdthis-threadsleep-for-and-gcc
#define _GLIBCXX_USE_NANOSLEEP 1

#include <gtest/gtest.h>

#include <iostream>
#include <array>
#include <algorithm>
#include <limits>
#include <numeric>
#include <chrono>
#include <thread>
#include <fstream>
#include <set>
#include <vector>

#include <networkit/auxiliary/ArrayTools.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/auxiliary/MissingMath.hpp>
#include <networkit/auxiliary/PrioQueue.hpp>
#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/auxiliary/SetIntersector.hpp>
#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/NumberParsing.hpp>
#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/BloomFilter.hpp>

namespace NetworKit {

class AuxGTest: public testing::Test{};

TEST_F(AuxGTest, testTimer) {
    Aux::LoggingTimer ltimer("LogginTimer", Aux::Log::LogLevel::trace); // only enabled if test is execut with raised loglevel

    Aux::StartedTimer timer;
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    const auto ms = timer.elapsedMilliseconds();
    const auto us = timer.elapsedMicroseconds();
    const auto ns = timer.elapsedNanoseconds();

    ASSERT_GE(ms, 5);
    ASSERT_GE(us, 5000);
    ASSERT_GE(ns, 5000000);

    timer.stop();

    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    const auto ms2 = timer.elapsedMilliseconds();
    const auto us2 = timer.elapsedMicroseconds();
    const auto ns2 = timer.elapsedNanoseconds();

    ASSERT_GE(ms2, 5);
    ASSERT_GE(us2, 5000);
    ASSERT_GE(ns2, 5000000);

    ASSERT_LE(ms2, us2 / 1000 + 1);
    ASSERT_LE(us2, ns2 / 1000 + 1);
}

TEST_F(AuxGTest, produceRandomIntegers) {
    Aux::Random::setSeed(1, false);
#ifndef NETWORKIT_RELEASE_LOGGING
    int64_t l = 0; 	// lower bound
    int64_t u = 100;	// upper bound

    for (int i = 0; i < 100; ++i) {
        TRACE(Aux::Random::integer(l, u));
    }
#else
    for (int i = 0; i < 100; ++i) {
        TRACE(Aux::Random::integer(1, 100));
    }
#endif
}

TEST_F(AuxGTest, testRandomInteger) {
    Aux::Random::setSeed(1, false);
    int64_t l = 0; 	// lower bound
    int64_t u = 10;	// upper bound
    std::vector<int64_t> rVector;
    int n = 1000;
    for (int i = 0; i < n; ++i) {
        int64_t r = Aux::Random::integer(l, u);
        assert(l <= r);
        assert(r <= u);
        rVector.push_back(r);
    }

    int64_t minR = *(min_element(rVector.begin(), rVector.end()));
    int64_t maxR = *(max_element(rVector.begin(), rVector.end()));

    EXPECT_EQ(minR, l);
    EXPECT_EQ(maxR, u);

    double sum = std::accumulate(rVector.begin(), rVector.end(), uint64_t{0});
    double avg = sum / n;


    DEBUG("avg rand integer: ", avg);
    EXPECT_LE(avg, 6.0);
    EXPECT_GE(avg, 4.0);
}

TEST_F(AuxGTest, testRandomReal) {
    Aux::Random::setSeed(1, false);
    std::vector<double> rVector;
    int n = 1000;
    for (int i = 0; i < n; ++i) {
        double r = Aux::Random::real();
        assert(0.0 <= r);
        assert(r < 1.0);
        rVector.push_back(r);
    }

    double sum = std::accumulate(rVector.begin(), rVector.end(), 0.0);
    double avg = sum / n;


    DEBUG("avg rand probability: ", avg);
    EXPECT_LE(avg, 0.6);
    EXPECT_GE(avg, 0.4);
}

TEST_F(AuxGTest, testRandomProbability) {
    Aux::Random::setSeed(1, false);
    std::vector<double> rVector;
    int n = 1000;
    for (int i = 0; i < n; ++i) {
        double r = Aux::Random::probability();
        assert(0.0 <= r);
        assert(r <= 1.0);
        rVector.push_back(r);
    }

    double sum = std::accumulate(rVector.begin(), rVector.end(), 0.0);
    double avg = sum / n;


    DEBUG("avg rand probability: ", avg);
    EXPECT_LE(avg, 0.6);
    EXPECT_GE(avg, 0.4);
}


TEST_F(AuxGTest, testBinomial) {
    EXPECT_EQ(1, Aux::MissingMath::binomial(7,0));
    EXPECT_EQ(7, Aux::MissingMath::binomial(7,1));
    EXPECT_EQ(21, Aux::MissingMath::binomial(7,2));
    EXPECT_EQ(35, Aux::MissingMath::binomial(7,3));
    EXPECT_EQ(35, Aux::MissingMath::binomial(7,4));
    EXPECT_EQ(21, Aux::MissingMath::binomial(7,5));
    EXPECT_EQ(7, Aux::MissingMath::binomial(7,6));
    EXPECT_EQ(1, Aux::MissingMath::binomial(7,7));


}


TEST_F(AuxGTest, benchmarkBinomial) {
    Aux::Timer timer;
    INFO("starting calculation");
    timer.start();

    int64_t n = 500;
    for (int64_t k = 0; k < n; ++k) {
        Aux::MissingMath::binomial(n, k);
    }

    timer.stop();

    INFO("calculation finished after ", timer.elapsedTag());

}

TEST_F(AuxGTest, testPriorityQueue) {
    using ElemType = std::pair<double, uint64_t>;

    // fill vector with keys
    std::vector<double> vec = {0.5, 3.5, 4.5, 2.5, 0.75, 1.5, 8.5, 3.25, 4.75, 5.0, 11.5, 0.25};

    // construct pq from vector
    Aux::PrioQueue<double, uint64_t> pq(vec);
    EXPECT_EQ(pq.size(), vec.size());
    EXPECT_FALSE(pq.empty());

    double topKey = pq.peekMin(0).first;
    pq.forElements([&](double curKey, uint64_t){
        EXPECT_TRUE(curKey >= topKey);
        topKey = curKey;
    });

    for (uint64_t elem = 0; elem < vec.size(); ++elem)
        EXPECT_TRUE(pq.contains(elem));
    EXPECT_FALSE(pq.contains(vec.size()));

    ElemType elem = pq.extractMin();
    EXPECT_EQ(0.25, elem.first);
    EXPECT_EQ(11u, elem.second);
    EXPECT_EQ(pq.size(), vec.size() - 1);
    EXPECT_FALSE(pq.empty());
    EXPECT_FALSE(pq.contains(11u));

    elem = pq.extractMin();
    EXPECT_EQ(0.5, elem.first);
    EXPECT_EQ(0u, elem.second);
    EXPECT_EQ(pq.size(), vec.size() - 2);
    EXPECT_FALSE(pq.empty());
    EXPECT_FALSE(pq.contains(0u));

    elem = pq.extractMin();
    EXPECT_EQ(0.75, elem.first);
    EXPECT_EQ(4u, elem.second);
    EXPECT_EQ(pq.size(), vec.size() - 3);
    EXPECT_FALSE(pq.empty());
    EXPECT_FALSE(pq.contains(4u));

    elem = pq.extractMin();
    EXPECT_EQ(1.5, elem.first);
    EXPECT_EQ(5u, elem.second);
    EXPECT_EQ(pq.size(), vec.size() - 4);
    EXPECT_FALSE(pq.empty());
    EXPECT_FALSE(pq.contains(5u));

    elem = pq.extractMin();
    EXPECT_EQ(2.5, elem.first);
    EXPECT_EQ(3u, elem.second);
    EXPECT_EQ(pq.size(), vec.size() - 5);
    EXPECT_FALSE(pq.empty());
    EXPECT_FALSE(pq.contains(3u));

    do {
        pq.extractMin();
    } while (pq.size() > 0);
    EXPECT_EQ(pq.size(), 0);
    EXPECT_TRUE(pq.empty());
    for (uint64_t i = 0; i < vec.size(); ++i)
        EXPECT_FALSE(pq.contains(i));
}

TEST_F(AuxGTest, testBucketPQ) {
    Aux::Random::setSeed(42, false);
    const int64_t maxKey = 100;
    const count capacity = 100;
    Aux::BucketPQ prioQ(capacity, 0, maxKey);

    EXPECT_TRUE(prioQ.empty());

    do {
        prioQ.insert(Aux::Random::integer(maxKey), prioQ.size());
        EXPECT_TRUE(prioQ.contains(prioQ.size() - 1));
        EXPECT_FALSE(prioQ.contains(prioQ.size()));
    } while (prioQ.size() < capacity);

    EXPECT_FALSE(prioQ.empty());
    for (index i = 0; i < prioQ.size(); ++i)
        EXPECT_TRUE(prioQ.contains(i));

    int64_t topKey = prioQ.extractMin().first;
    do {
        EXPECT_GE(prioQ.getMin().first, topKey);
        topKey = prioQ.extractMin().first;
    } while (!prioQ.empty());

    EXPECT_TRUE(prioQ.empty());
    for (index i = 0; i < prioQ.size(); ++i)
        EXPECT_FALSE(prioQ.contains(i));
}

TEST_F(AuxGTest, testBucketPQUpdateRemove) {
    Aux::Random::setSeed(42, false);
    const int64_t maxKey = 100;
    const count capacity = 100;
    Aux::BucketPQ prioQ(capacity, 0, maxKey);

    EXPECT_TRUE(prioQ.empty());

    do {
        prioQ.insert(Aux::Random::integer(maxKey), prioQ.size());
    } while (prioQ.size() < capacity);

    EXPECT_FALSE(prioQ.empty());

    for (count i = 0; i < capacity; ++i) {
        EXPECT_TRUE(prioQ.contains(i));
        const auto newKey = Aux::Random::integer(maxKey);
        prioQ.changeKey(newKey, i);
        EXPECT_EQ(prioQ.getKey(i), newKey);
        EXPECT_TRUE(prioQ.contains(i));
    }

    for (count i = 0; i < capacity; ++i) {
        EXPECT_FALSE(prioQ.empty());
        prioQ.remove(i);
        EXPECT_FALSE(prioQ.contains(i));
    }

    EXPECT_TRUE(prioQ.empty());
}

TEST_F(AuxGTest, testLogging) {
    std::string cl = Aux::Log::getLogLevel();
    //FATAL("FATAL ERROR");
    //ERROR("This may be here");
    Aux::Log::setLogLevel("ERROR");
    EXPECT_STREQ("ERROR",Aux::Log::getLogLevel().c_str());
    //FATAL("fatal error",3,2.f);
    //ERROR("normal error",3,2.f);
    //WARN("just a warning",3,2.f);
    //INFO("for the sake of information", 3, 2.f);
    //DEBUG("important debug outputs", 3, 2.f);
    //TRACE("trace", 3, 2.f);
    Aux::Log::setLogLevel("INFO");
    EXPECT_STREQ("INFO",Aux::Log::getLogLevel().c_str());
    //FATAL("fatal error", 3, 2.f);
    //ERROR("normal error", 3, 2.f);
    //WARN("just a warning", 3, 2.f);
    //INFO("for the sake of information", 3, 2.f);
    //DEBUG("important debug outputs", 3, 2.f);
    //TRACE("trace", 3, 2.f);
    Aux::Log::setLogLevel("TRACE");
    EXPECT_STREQ("TRACE",Aux::Log::getLogLevel().c_str());
    //FATAL("fatal error", 3, 2.f);
    //ERROR("normal error", 3, 2.f);
    //WARN("just a warning", 3, 2.f);
    //INFO("for the sake of information", 3, 2.f);
    //DEBUG("important debug outputs", 3, 2.f);
    //TRACE("trace", 3, 2.f);
    Aux::Log::setLogLevel("WARN");
    EXPECT_STREQ("WARN",Aux::Log::getLogLevel().c_str());
    Aux::Log::setLogLevel("FATAL");
    EXPECT_STREQ("FATAL",Aux::Log::getLogLevel().c_str());
    Aux::Log::setLogLevel("DEBUG");
    EXPECT_STREQ("DEBUG",Aux::Log::getLogLevel().c_str());
    Aux::Log::setLogLevel(cl);
}

TEST_F(AuxGTest, testFormatting) {
    using Aux::toStringF;
    EXPECT_EQ("",toStringF(""));
    EXPECT_EQ("",toStringF("%s", ""));
    EXPECT_EQ("%",toStringF("%%", ""));
    EXPECT_EQ("aaabbbaaa",toStringF("aaa%saaa", "bbb"));
    EXPECT_THROW(toStringF("%"), std::invalid_argument);
    EXPECT_THROW(toStringF("%i"), std::invalid_argument);
    EXPECT_THROW(toStringF("%s"), std::invalid_argument);
}

TEST_F(AuxGTest, testRandomChoice) {
    std::vector<uint64_t> data = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    for (uint64_t i = 0; i < 1000; ++i) {
        EXPECT_TRUE(std::find(data.begin(), data.end(), Aux::Random::choice(data)) != data.end());
    }
}


TEST_F(AuxGTest, testRandomWeightedChoice) {
    std::vector<std::pair<uint64_t, double> > data = {{0, 1.0}, {1, 0.0}};
    for (uint64_t i = 0; i < 100; ++i) {
        auto element = Aux::Random::weightedChoice(data);
        EXPECT_EQ(0u, element);
    }
}

TEST_F(AuxGTest, testRandomIndex) {
    Aux::Random::setSeed(1, false);

    for (unsigned i = 0; i < 10; i++) {
        EXPECT_EQ(0u, Aux::Random::index(1));
    }

    for (unsigned i = 0; i < 100; i++) {
        auto tmp = Aux::Random::index(10);
        EXPECT_LE(tmp, 9u);
        EXPECT_GE(tmp, 0u);
    }
}

TEST_F(AuxGTest, testSplit) {
    using Vec = std::vector<std::string>;

    EXPECT_EQ(Vec{}, Aux::StringTools::split(""));
    EXPECT_EQ(Vec{""}, Aux::StringTools::split(" "));

    {
        auto expected = Vec{"", ""};
        EXPECT_EQ(expected, Aux::StringTools::split("  "));
    }
    {
        auto expected = Vec{"", "a"};
        EXPECT_EQ(expected, Aux::StringTools::split(" a"));
    }
    {
        auto expected = Vec{"a"};
        EXPECT_EQ(expected, Aux::StringTools::split("a "));
    }
    {
        auto expected = Vec{"a"};
        EXPECT_EQ(expected, Aux::StringTools::split("a"));
    }
    {
        auto expected = Vec{"a", "b"};
        EXPECT_EQ(expected, Aux::StringTools::split("a b"));
    }
    {
        auto expected = Vec{"", "a", "b"};
        EXPECT_EQ(expected, Aux::StringTools::split(" a b "));
    }
    {
        auto expected = Vec{"abc", "def", "ghi"};
        EXPECT_EQ(expected, Aux::StringTools::split("abc def ghi"));
    }
    {
        auto expected = Vec{"abc", "def", "ghi"};
        EXPECT_EQ(expected, Aux::StringTools::split("abc def ghi "));
    }
}


TEST_F(AuxGTest, testSetIntersector) {
    uint64_t n = 50;
    std::vector<uint64_t> A = {0, 5, 7, 3, 11, 22, 45, 17, 23, 11, 45};
    std::vector<uint64_t> B = {27, 11, 13, 17, 44, 1, 2, 3, 9, 44, 11, 20};

    Aux::SetIntersector<uint64_t> si(n);
    std::set<uint64_t> intersection = si.intersect(A, B);

    std::set<uint64_t> expectedResult = {3, 11, 17};
    EXPECT_EQ(expectedResult, intersection);
}

TEST_F(AuxGTest, testEnforce) {
    EXPECT_THROW(Aux::enforce(false), std::runtime_error);
    EXPECT_NO_THROW(Aux::enforce(true));

    EXPECT_THROW(Aux::enforce(false, "foo"), std::runtime_error);
    EXPECT_NO_THROW(Aux::enforce(true, "bar"));

    EXPECT_THROW(Aux::enforce<std::logic_error>(false, "foo"), std::logic_error);
    EXPECT_NO_THROW(Aux::enforce<std::logic_error>(true, "foo"));

    std::string msg = "some message in a std::string";

    EXPECT_THROW(Aux::enforce(false, msg), std::runtime_error);
    EXPECT_NO_THROW(Aux::enforce(true, msg));

    EXPECT_THROW(Aux::enforce<std::logic_error>(false, msg), std::logic_error);
    EXPECT_NO_THROW(Aux::enforce<std::logic_error>(true, msg));
}

TEST_F(AuxGTest, testEnforceOpened) {
    using Aux::enforceOpened;
    std::ifstream stream;
    EXPECT_THROW(enforceOpened(stream), std::runtime_error);
}

TEST_F(AuxGTest, testNumberParsingInteger) {
    const std::string str = "0 00 1 123 001 1200 12345678    ";
    std::vector<unsigned> expectedValues = {0, 0, 1, 123, 1, 1200, 12345678};
    auto it = str.begin();
    auto end = str.end();
    std::size_t i = 0;
    while(it != end) {
        unsigned result;
        std::tie(result, it) = Aux::Parsing::strTo<unsigned>(it, end);
        EXPECT_EQ(expectedValues[i], result);
        ++i;
    }

    EXPECT_EQ(it, end);
}

TEST_F(AuxGTest, testNumberParsingSignedInteger) {
    const std::string str = "-0 -00 -1 -123 -001 -1200 -12345678    ";
    std::vector<int> expectedValues = {0, 0, -1, -123, -1, -1200, -12345678};
    auto it = str.begin();
    auto end = str.end();
    std::size_t i = 0;
    while(it != end) {
        int result;
        std::tie(result, it) = Aux::Parsing::strTo<int>(it, end);
        EXPECT_EQ(expectedValues[i], result);
        ++i;
    }

    EXPECT_EQ(it, end);
}

TEST_F(AuxGTest, testOverflowCatching) {
    const std::string str = "1000";
    EXPECT_THROW(
            (Aux::Parsing::strTo<uint8_t, std::string::const_iterator, Aux::Checkers::Enforcer>(
                str.begin(), str.end())),
            std::runtime_error
            );
    EXPECT_THROW(
            (Aux::Parsing::strTo<int8_t, std::string::const_iterator, Aux::Checkers::Enforcer>(
                str.begin(), str.end())),
            std::runtime_error
            );
}

TEST_F(AuxGTest, testNumberParsingBasicReal) {
    const std::string str =
        "0 00 1 123 001 1200 12345678    "
        "0.00000 -0000.000 -0000.000e-100"
        ;
    std::vector<double> expectedValues = {
        0, 0, 1, 123, 1, 1200, 12345678,
        0, 0, 0
    };
    auto it = str.begin();
    auto end = str.end();
    std::size_t i = 0;
    while(it != end) {
        double result;
        std::tie(result, it) = Aux::Parsing::strTo<double>(it, end);
        EXPECT_EQ(expectedValues[i], result);
        ++i;
    }

    EXPECT_EQ(it, end);
}


TEST_F(AuxGTest, testNumberParsingAdvancedReal) {
    auto helper = [](const std::string& str, double expected) {
        auto result = std::get<0>(Aux::Parsing::strTo<double>(str.begin(), str.end()));
        EXPECT_DOUBLE_EQ(result, expected);
    };
#define TEST_CASE_REAL(number) helper(#number, number)
    TEST_CASE_REAL(0);
    TEST_CASE_REAL(1);
    TEST_CASE_REAL(8.76464e+23);
    TEST_CASE_REAL(5.56628e+12);
    TEST_CASE_REAL(9.563574701896000000e+23);
    TEST_CASE_REAL(9.563574701896270940e+23);
    TEST_CASE_REAL(-669585235219883571019776.);
    TEST_CASE_REAL(9.563574701896270946e+23);
    TEST_CASE_REAL(141265720771594800000000.);
    TEST_CASE_REAL(141265720771594850000000.);
    TEST_CASE_REAL(1.4126572077159485e+23);
    TEST_CASE_REAL(1.4126572077159480e+23);
    TEST_CASE_REAL(-784008854951861199831040.);
    TEST_CASE_REAL(-78400885495186119983104.0);
#undef TEST_CASE_REAL
}

TEST_F(AuxGTest, testBloomFilter) {
    Aux::Random::setSeed(1, false);
    Aux::BloomFilter bf(5);
    uint64_t size = 750;
    std::set<uint64_t> randomKeys;

    while (randomKeys.size() < size) {
        uint64_t key = Aux::Random::integer();
        bf.insert(key);
        randomKeys.insert(key);
    }

    for (auto key: randomKeys) {
        EXPECT_TRUE(bf.isMember(key));
    }

    for (uint64_t i = 0; i < size; ++i) {
        uint64_t key = Aux::Random::integer();
        if (randomKeys.count(key) == 0) {
            if (bf.isMember(key)) {
                WARN("Bloom filter: Maybe only false positive, maybe indication of bug!");
            }
        }
    }
}

template <class Perm>
void testPermutation(Perm perm) {
    auto &urng = Aux::Random::getURNG();
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), urng);

    std::vector<index> seq;
    seq.resize(perm.size());
    std::for_each(seq.begin(), seq.end(), [&seq](auto &elem) { elem = Aux::Random::integer(seq.size()); });
    auto seqCopy = seq;

    Aux::ArrayTools::applyPermutation(seq.begin(), seq.end(), perm.begin());
    for (size_t i = 0; i < seq.size(); ++i)
        EXPECT_EQ(seqCopy[perm[i]], seq[i]);
}

TEST_F(AuxGTest, testApplyPermutation) {
    Aux::Random::setSeed(1, true);
    constexpr int n = 100;

    testPermutation(std::array<int8_t, n>{});
    testPermutation(std::array<int8_t, 128>{});
    testPermutation(std::array<uint8_t, n>{});
    testPermutation(std::array<uint8_t, 256>{});
    testPermutation(std::array<int16_t, n>{});
    testPermutation(std::array<int16_t, static_cast<size_t>(std::numeric_limits<int16_t>::max()) + 1>{});
    testPermutation(std::array<uint16_t, n>{});
    testPermutation(std::array<uint16_t, static_cast<size_t>(std::numeric_limits<uint16_t>::max()) + 1>{});
    testPermutation(std::array<uint32_t, n>{});
    testPermutation(std::array<int, n>{});
    testPermutation(std::array<int64_t, n>{});
    testPermutation(std::vector<node>{n});
    testPermutation(std::vector<int>{n});
}
} // ! namespace NetworKit
