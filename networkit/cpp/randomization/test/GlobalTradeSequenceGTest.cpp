/*
 * GlobalTradeSequenceGTest.cpp
 *
 *  Created on: 19.06.2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#include <map>
#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>

#include "../GlobalTradeSequence.hpp"

namespace NetworKit {

template <typename T>
class GlobalTradeSequenceGTest : public ::testing::Test {};

using HashMapTypes = ::testing::Types<CurveballDetails::FixedLinearCongruentialMap<uint64_t>,
                                      CurveballDetails::LinearCongruentialMap<uint64_t>>;
TYPED_TEST_SUITE(GlobalTradeSequenceGTest, HashMapTypes, );

TYPED_TEST(GlobalTradeSequenceGTest, testInvert) {
    using Hash = TypeParam;
    using value_type = typename Hash::value_type;

    Aux::Random::setSeed(1, false);

    for (int repeat = 0; repeat < 10; repeat++) {
        for (value_type n : {10, 123, 234, 3456, 123456}) {
            Hash hash(n, Aux::Random::getURNG());

            for (value_type i = 0; i < n; i++) {
                ASSERT_EQ(hash.invert(hash(i)), i) << "i=" << i;
            }
        }
    }
}

TYPED_TEST(GlobalTradeSequenceGTest, testInjectivity) {
    using Hash = TypeParam;
    using value_type = typename Hash::value_type;

    Aux::Random::setSeed(1, false);

    for (int repeat = 0; repeat < 10; repeat++) {
        for (value_type n : {10, 123, 234, 3456, 123456}) {
            std::map<value_type, value_type> image;
            Hash hash(n, Aux::Random::getURNG());

            for (value_type i = 0; i < n; i++) {
                const auto hashed = hash(i);

                const auto it = image.find(hashed);
                ASSERT_TRUE(it == image.end())
                    << "i=" << i << ", hashed=" << hashed << ", earlier=" << it->second;

                image[hashed] = i;
            }
        }
    }
}

} // namespace NetworKit
