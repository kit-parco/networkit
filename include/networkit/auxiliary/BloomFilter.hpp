// no-networkit-format
/*
 * BloomFilter.h
 *
 *  Created on: 08.08.2015
 *      Author: Henning
 */

#ifndef NETWORKIT_AUXILIARY_BLOOM_FILTER_HPP_
#define NETWORKIT_AUXILIARY_BLOOM_FILTER_HPP_

#include <vector>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/Globals.hpp>

namespace Aux {

using index = NetworKit::index;
using count = NetworKit::count;

/**
 * Bloom filter for fast set membership queries for type index with one-sided error (false positives).
 * Uses one hash function instead of many but different random salts instead.
 */
class BloomFilter {
private:
    count numHashes;
    count size;
    std::vector<std::vector<char> > membership;
    std::vector<index> salts;

    /**
     * Hashes @a key with the salt at index @a hfunc (to mimic a different hash function).
     * @param[in] key Key to be hashed.
     * @param[in] hfunc In principle index of a hash function, here index of a salt.
     * @return Hash value of salted key.
     */
    index hash(index key, index hfunc) const;

public:
    /**
     * Constructor. Space complexity: @a numHashes * @a size bytes.
     * @param[in] numHashes Number of different hash functions that should be used.
     * Actually, the implementation uses @a numHashes different random salts instead and one
     * hash function only.
     * @param[in] size Size of each hash array.
     */
    BloomFilter(count numHashes, count size = 6291469);

    /**
     * Destructor.
     */
    virtual ~BloomFilter() = default;

    /**
     * Inserts @a key into Bloom filter.
     * @param[in] key Key to be inserted.
     */
    void insert(index key);

    /**
     * @param[in] key Key to be queried for set membership.
     * @return True if @a is member, false otherwise.
     */
    bool isMember(index key) const;
};

}

#endif // NETWORKIT_AUXILIARY_BLOOM_FILTER_HPP_
