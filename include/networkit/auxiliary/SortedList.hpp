// no-networkit-format
/*
 * SortedList.h
 *
 *  Created on: 21.09.2018
 *      Author: Eugenio Angriman
 */

#ifndef NETWORKIT_AUXILIARY_SORTED_LIST_HPP_
#define NETWORKIT_AUXILIARY_SORTED_LIST_HPP_

#include <algorithm>
#include <vector>

#include <networkit/auxiliary/Log.hpp>

namespace Aux {
/*
 * Keeps a sorted list of pairs with at most k elements.
 * If more than k elements are inserted, the elements with smallest value are
 * removed. The list is implemented on top of vector, thus the insert operation
 * takes O(k) time.
 *
 * Warning: this sorted list was designed for the Kadabra algorithm; we assume
 * that all the newly inserted elements have a value that is greater or equal
 * to 0 if not already present in the list, or greater or equal their previous
 * value.
 * TODO: generalize this data structure.
 */
class SortedList {
private:
    std::vector<std::pair<uint64_t, double>> elements;
    std::vector<uint64_t> position;
    uint64_t virtualSize;
    const uint64_t size;
    const uint64_t maxKey;

public:
    /**
     * Creates a SortedList of size @a size accepting key values from 0 to
     * @a maxKey.
     */
    SortedList(uint64_t size, uint64_t maxKey);

    /**
     * Insert a key-value element.
     */
    void insert(uint64_t newElement, double newValue);

    /**
     * Returns the value at position @a i in the ranking.
     */
    double getValue(const uint64_t i) const { return elements[i].second; }

    /**
     * Returns the element at position @a i in the ranking.
     */
    uint64_t getElement(const uint64_t i) const { return elements[i].first; }

    /**
     * Returns the number of elements stored in the list.
     */
    uint64_t getSize() const { return virtualSize; }

    /**
     * Removes all the elements in the list.
     */
    void clear();
};

inline SortedList::SortedList(const uint64_t size, const uint64_t maxKey)
    : size(size), maxKey(maxKey) {
    if (maxKey < size) {
        throw std::runtime_error("maxKey must be bigger than the size.");
    }
    clear();
}

inline void SortedList::clear() {
    elements.resize(size);
    position.resize(maxKey);
    uint64_t i;
    for (i = 0; i < size; ++i) {
        elements[i] = std::make_pair(i, 0.);
        position[i] = i;
    }
    for (i = size; i < maxKey; ++i) {
        position[i] = i;
    }
    virtualSize = 0;
}

inline void SortedList::insert(const uint64_t newElement, const double newValue) {
    uint64_t ub =
        std::upper_bound(
            elements.begin(), elements.begin() + virtualSize, newValue,
            [&](const double val, const std::pair<uint64_t, double> pair) {
                return val > pair.second;
            }) -
        elements.begin();

    uint64_t oldPos;
    // We assume that if the same key is inserted again, its value will be
    // greater than before (i.e., oldPos >= ub).
    if (ub < size) {
        oldPos = std::min(position[newElement], size - 1);
        if (position[newElement] < size) {
            assert(elements[ub].second <= newValue);
        }
        if (virtualSize < size && oldPos >= virtualSize) {
            ++virtualSize;
        }

        if (ub < oldPos) {
            std::rotate(elements.begin() + ub, elements.begin() + oldPos,
                        elements.begin() + oldPos + 1);

            if (oldPos == size - 1) {
                ++position[elements[ub].first];
            }

            elements[ub] = std::make_pair(newElement, newValue);
            for (auto it = elements.begin() + ub + 1;
                 it < elements.begin() + oldPos + 1; ++it) {
                ++position[(*it).first];
            }
        } else {
            elements[ub].second = newValue;
        }
    }
}
} // namespace Aux

#endif // NETWORKIT_AUXILIARY_SORTED_LIST_HPP_
