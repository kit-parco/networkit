// no-networkit-format
/*
 * BucketPQ.cpp
 *
 *  Created on: 03.03.2017
 *      Author: Henning
 */

#include <networkit/auxiliary/BucketPQ.hpp>

namespace Aux {

Bucket BucketPQ::dummyBucket = {};
const Bucket::iterator BucketPQ::invalidPtr = BucketPQ::dummyBucket.end();

BucketPQ::BucketPQ(const std::vector<int64_t>& keys, int64_t minAdmissibleKey, int64_t maxAdmissibleKey):
        minAdmissibleKey(minAdmissibleKey), maxAdmissibleKey(maxAdmissibleKey)
{
    construct(keys.size());

    // insert key-value pairs
    for (index i = 0; i < keys.size(); ++i) {
        if (keys[i] != none) {
            insert(keys[i], i);
        }
    }
}

BucketPQ::BucketPQ(uint64_t capacity, int64_t minAdmissibleKey, int64_t maxAdmissibleKey):
        minAdmissibleKey(minAdmissibleKey), maxAdmissibleKey(maxAdmissibleKey)
{
    construct(capacity);
}

void BucketPQ::construct(uint64_t capacity) {
    // check range of keys
    if (minAdmissibleKey > maxAdmissibleKey) {
        throw std::invalid_argument("minAdmissibleKey cannot be larger than maxAdmissibleKey");
    }

    // init
    buckets.resize(maxAdmissibleKey-minAdmissibleKey+1);
    nodePtr.resize(capacity);
    myBucket.resize(capacity);
    currentMinKey = std::numeric_limits<int64_t>::max();
    currentMaxKey = std::numeric_limits<int64_t>::min();
    numElems = 0;

    offset = -minAdmissibleKey;
}

void BucketPQ::insert(int64_t key, index value) {
    assert(minAdmissibleKey <= key && key <= maxAdmissibleKey);
    assert(value < nodePtr.size());

    buckets[key+offset].push_front(value);
    nodePtr[value] = OptionalIterator{true, buckets[key+offset].begin()};
    myBucket[value] = key+offset;
    ++numElems;

    // bookkeeping
    if (key < currentMinKey) {
        currentMinKey = key;
    }
    if (key > currentMaxKey) {
        currentMaxKey = key;
    }
}

bool BucketPQ::contains(const index &value) const {
    return value < nodePtr.size() && nodePtr[value].valid;
}

void BucketPQ::remove(const index& value) {
    assert(value < nodePtr.size());

    if (myBucket[value] != none) {
        // remove from appropriate bucket
        index key = myBucket[value];
        buckets[key].erase(nodePtr[value].iter);
        nodePtr[value].reset();
        myBucket[value] = none;
        --numElems;

        if (size() == 0) {
            // empty pq: reinit the current min/max pointers
            currentMinKey = std::numeric_limits<int64_t>::max();
            currentMaxKey = std::numeric_limits<int64_t>::min();
        }
        else {
            // adjust max pointer if necessary
            while (buckets[currentMaxKey+offset].empty() && currentMaxKey > currentMinKey) {
                --currentMaxKey;
            }

            // adjust min pointer if necessary
            while (buckets[currentMinKey+offset].empty() && currentMinKey < currentMaxKey) {
                ++currentMinKey;
            }
        }
    }
}

std::pair<int64_t, index> BucketPQ::extractMin() {
    if (empty())
        return {none, none};

    index result = buckets[currentMinKey+offset].front();

    // store currentMinKey because remove(result) will change it
    int64_t oldMinKey = currentMinKey;
    remove(result);
    return {oldMinKey, result};
}

std::pair<int64_t, index> BucketPQ::getMin() {
    if (empty())
        return {none, none};
    else
        return {currentMinKey, buckets[currentMinKey + offset].front()};
}

void BucketPQ::changeKey(int64_t newKey, index value) {
    remove(value);
    insert(newKey, value);
}

uint64_t BucketPQ::size() const {
    return numElems;
}

bool BucketPQ::empty() const noexcept {
    return numElems == 0;
}

int64_t BucketPQ::getKey(const index& val) {
    int64_t key = static_cast<int64_t>(myBucket[val]) - offset;
    return key;
}
} // namespace
