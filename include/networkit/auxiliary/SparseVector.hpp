/*
 * SparseVector.hpp
 *
 * Created: 2019-10-15
 * Author: Armin Wiebigke
 */
#ifndef NETWORKIT_AUXILIARY_SPARSE_VECTOR_HPP_
#define NETWORKIT_AUXILIARY_SPARSE_VECTOR_HPP_

#include <algorithm>
#include <cassert>
#include <vector>

#include <networkit/Globals.hpp>

namespace NetworKit {

/**
 * A vector that imitates a map with unsigned integer keys.
 * This class has faster access than a map, but needs space linear in the maximum key value.
 */
template <typename T>
class SparseVector {
public:
    SparseVector();

    /**
     * Construct an empty vector. Empty values are created using the default constructor.
     * @param size upper bound for the maximum usable index
     */
    explicit SparseVector(index size);

    /**
     * Construct an empty vector.
     * @param size upper bound for the maximum usable index
     * @param emptyValue value used for empty entries
     */
    SparseVector(index size, T emptyValue);

    /**
     * Resize the vector so that indexes up to size-1 can be used.
     */
    void setUpperBound(index size);

    /**	 *
     * @return the upper bound of the indexes that can be used
     */
    index upperBound() const;

    /**
     * @return the number of inserted elements.
     */
    count size() const;

    /**
     * Insert a value at a given position.
     * @param i index where the value is inserted
     * @param value
     */
    void insert(index i, T value);

    /**
     * Access operator. Before accessing an element, insert it by using the insert() method.
     */
    T &operator[](index i);

    /**
     * Const access operator. Before accessing an element, insert it by using the insert() method.
     */
    const T &operator[](index i) const;

    /**
     * Returns true iff an element was previously inserted at the given index.
     * @param idx
     */
    bool indexIsUsed(index idx);

    /**
     * Inserts value at position i, or replaces the value if previously inserted
     * @param i
     * @param value
     */
    void insertOrSet(index i, T value);

    /**
     * Remove all indexes for which the value is set to the emptyValue.
     */
    void removeUnusedIndexes() {
        auto new_end = std::remove_if(usedIndexes.begin(), usedIndexes.end(),
                                      [&](index i) { return data[i] == emptyValue; });
        usedIndexes.erase(new_end, usedIndexes.end());
    }

    /**
     * Reset all values to the default value, so it is "empty". The upper bound is not changed.
     */
    void reset();

    /**
     * Clear the vector, setting the upper bound of usable indexes to 0.
     */
    void clear();

    /**
     * Reallocate the datastructure if size exceeds current upper bound
     * This is different from setUpperBound() since we want to make sure both usedIndexes and data
     * are allocated on the socket of the calling thread
     * @param size
     * @param emptyValue new emptyValue
     */
    void resize(size_t size, T emptyValue);

    /**
     * Applies the given lambda to each inserted index and associated value
     * @tparam ElementHandler
     * @param lambda
     */
    template <typename ElementHandler>
    void forElements(ElementHandler &&lambda) {
        for (index i : usedIndexes) {
            lambda(i, data[i]);
        }
    }

private:
    std::vector<T> data;
    std::vector<index> usedIndexes;
    T emptyValue;
};

template <typename T>
SparseVector<T>::SparseVector() : emptyValue(T{}) {}

template <typename T>
SparseVector<T>::SparseVector(count size) : SparseVector(size, T{}) {}

template <typename T>
SparseVector<T>::SparseVector(count size, T emptyValue)
    : data(size, emptyValue), emptyValue(emptyValue) {}

template <typename T>
void SparseVector<T>::reset() {
    for (index i : usedIndexes) {
        data[i] = emptyValue;
    }
    usedIndexes.clear();
}

template <typename T>
void SparseVector<T>::insert(index i, T value) {
    assert(data[i] == emptyValue);
    usedIndexes.push_back(i);
    data[i] = std::move(value);
}

template <typename T>
T &SparseVector<T>::operator[](index i) {
    return data[i];
}

template <typename T>
const T &NetworKit::SparseVector<T>::operator[](NetworKit::index i) const {
    assert(i < data.size());
    return data[i];
}

template <typename T>
void SparseVector<T>::setUpperBound(index size) {
    data.resize(size, emptyValue);
}

template <typename T>
index SparseVector<T>::upperBound() const {
    return data.size();
}

template <typename T>
count SparseVector<T>::size() const {
    return usedIndexes.size();
}

template <typename T>
void NetworKit::SparseVector<T>::clear() {
    usedIndexes.clear();
    usedIndexes.shrink_to_fit();
    data.clear();
    data.shrink_to_fit();
}

template <typename T>
void NetworKit::SparseVector<T>::resize(size_t size, T emptyValue) {
    if (size > upperBound()) {
        this->emptyValue = emptyValue;
        data = std::vector<T>(size, this->emptyValue);
        usedIndexes = std::vector<index>();
    }
}

template <typename T>
bool NetworKit::SparseVector<T>::indexIsUsed(index idx) {
    return data[idx] != emptyValue;
}

template <typename T>
void NetworKit::SparseVector<T>::insertOrSet(NetworKit::index i, T value) {
    if (!indexIsUsed(i)) {
        insert(i, value);
    } else {
        data[i] = value;
    }
}

} /* namespace NetworKit */

#endif // NETWORKIT_AUXILIARY_SPARSE_VECTOR_HPP_
