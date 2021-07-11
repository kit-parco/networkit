// no-networkit-format
/*
 * BucketPQ.hpp
 *
 *  Created on: 02.03.2017
 *      Author: Henning
 */

#ifndef NETWORKIT_AUXILIARY_BUCKET_PQ_HPP_
#define NETWORKIT_AUXILIARY_BUCKET_PQ_HPP_

#include <limits>
#include <list>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/PrioQueue.hpp>

namespace Aux {

using index = NetworKit::index;
using count = NetworKit::count;
using Bucket = std::list<index>;
constexpr int64_t none = std::numeric_limits<int64_t>::max();

/**
 * Addressable priority queue for values in the range [0,n) and
 * integer keys (= priorities) in the range [minPrio, maxPrio].
 * minPrio and maxPrio can be positive or negative, respectively with
 * the obvious constraint minPrio <= maxPrio.
 * Amortized constant running time for each operation.
 */
class BucketPQ: public PrioQueue<int64_t, index> {
private:
    std::vector<Bucket> buckets;			// the actual buckets
    static Bucket dummyBucket;
    static const Bucket::iterator invalidPtr;

    struct OptionalIterator {
        bool valid;
        Bucket::iterator iter;

        void reset() {
            valid = false;
            iter = invalidPtr;
        }
        OptionalIterator() { reset(); }
        OptionalIterator(bool valid, Bucket::iterator iter) : valid(valid), iter(iter) {}
    };

    std::vector<OptionalIterator> nodePtr; // keeps track of node positions
    std::vector<index> myBucket;			// keeps track of current bucket for each value
    int64_t currentMinKey;					// current min key
    int64_t currentMaxKey;					// current max key
    int64_t minAdmissibleKey;				// minimum admissible key
    int64_t maxAdmissibleKey;				// maximum admissible key
    count numElems;							// number of elements stored
    int64_t offset;							// offset from minAdmissibleKeys to 0

    /**
     * Constructor. Not to be used, only here for overriding.
     */
    BucketPQ(const std::vector<int64_t>&) {}

    /**
     * Constructor. Not to be used, only here for overriding.
     */
    BucketPQ(uint64_t) {}

    /**
     * Called from various constructors for initializing members.
     */
    void construct(uint64_t capacity);

public:

    /**
     * Builds priority queue from the vector @a keys, values are indices
     * of @a keys.
     * @param[in] keys Vector of keys
     * @param[in] minAdmissibleKey Minimum admissible key
     * @param[in] maxAdmissibleKey Maximum admissible key
     */
    BucketPQ(const std::vector<int64_t>& keys, int64_t minAdmissibleKey, int64_t maxAdmissibleKey);

    /**
    * Builds priority queue of the specified capacity @a capacity.
    */
    BucketPQ(uint64_t capacity, int64_t minAdmissibleKey, int64_t maxAdmissibleKey);

    /**
     * Default destructor
     */
    ~BucketPQ() override = default;

    /**
     * Inserts key-value pair (@key, @value).
     */
    void insert(int64_t key, index value) override;

    /**
     * Returns the element on top of the PrioQ.
     */
    std::pair<int64_t, index> getMin();

    /**
     * Removes the element with minimum key and returns the key-value pair.
     */
    std::pair<int64_t, index> extractMin() override;

    /**
     * Modifies entry with value @a value.
     * The entry is then set to @a newKey with the same value.
     * If the corresponding key is not present, the element will be inserted.
     */
    void changeKey(int64_t newKey, index value) override;

    /**
     * @return Number of elements in PQ.
     */
    uint64_t size() const override;

    /**
     * @return Whether or not the PQ is empty.
     */
    bool empty() const noexcept override;

    /**
     * @return Whether or not the PQ contains the given element.
     */
    bool contains(const index &value) const override;

    /**
     * Removes key-value pair given by value @a val.
     */
    void remove(const index &val) override;

    /**
     * @return key to given value @val.
     */
    virtual int64_t getKey(const index& val);

};

} /* namespace Aux */
#endif // NETWORKIT_AUXILIARY_BUCKET_PQ_HPP_
