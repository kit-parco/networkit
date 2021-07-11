// no-networkit-format
/*
 * BloomFilter.cpp
 *
 *  Created on: 08.08.2015
 *      Author: Henning
 */

#include <networkit/auxiliary/BloomFilter.hpp>
#include <networkit/auxiliary/Random.hpp>

namespace Aux {

BloomFilter::BloomFilter(count numHashes, count size): numHashes(numHashes), size(size),
    membership(numHashes), salts(numHashes)
{
    for (index i = 0; i < numHashes; ++i) {
        membership[i].resize(size, false);
        salts[i] = Aux::Random::integer();
    }
}

void BloomFilter::insert(index key) {
    // set hashed positions of key in each array to true
    for (index func = 0; func < numHashes; ++func) {
        membership[func][hash(key, func)] = true;
    }
}

bool BloomFilter::isMember(index key) const {
    bool result = true;

    // find out if all hashed positions of key are true
    for (index func = 0; (func < numHashes) && result; ++func) {
        result = result && membership[func][hash(key, func)];
    }
    return result;
}

index BloomFilter::hash(index key, index hfunc) const {
    index result = 0;
    std::hash<index> myhash;

    // instead of using different hash functions: salt key with different (but fixed) random salt
    result = myhash(key ^ salts[hfunc]);
    return (result % size);
}

}
