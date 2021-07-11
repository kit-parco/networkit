// no-networkit-format
/*
 * SetIntersector.h
 *
 *  Created on: 14.05.2014
 *      Author: Henning
 */

#ifndef NETWORKIT_AUXILIARY_SET_INTERSECTOR_HPP_
#define NETWORKIT_AUXILIARY_SET_INTERSECTOR_HPP_

#include <cstdint>
#include <vector>
#include <set>

namespace Aux {

/**
 * Provides set intersection for sets with entries from 0 to an upper bound (exclusive).
 * Preprocessing time: Linear time for object creation (creates a bit vector of size upperBound.
 * Query time: O(|A|+|B|) for two sets A and B.
 * Space complexity: upperBound + 64 bits.
 */
template<class T>
class SetIntersector {
public:
    /**
     * @param upperBound Exclusive upper bound for IDs of set members.
     */
    SetIntersector(T upperBound);

    /**
     * @return Intersection of sets provided in @a A and @a B.
     */
    std::set<T> intersect(const std::vector<T>& A, const std::vector<T>& B);

private:
    std::vector<bool> bv;
    uint64_t n;
};

} // namespace Aux


template<class T>
inline Aux::SetIntersector<T>::SetIntersector(T upperBound): n(upperBound) {
    bv.resize(n, false);

}

template<class T>
inline std::set<T> Aux::SetIntersector<T>::intersect(const std::vector<T>& A,
        const std::vector<T>& B)
{
    const std::vector<T>& smaller = (A.size() <= B.size()) ? A : B;
    const std::vector<T>& larger = (A.size() <= B.size()) ? B : A;

    for (auto entry: smaller) {
        bv[entry] = true;
    }

    std::set<T> result;
    for (auto entry: larger) {
        if (bv[entry]) {
            result.insert(entry);
        }
    }

    for (auto entry: smaller) {
        bv[entry] = false;
    }

    return result;
}

#endif // NETWORKIT_AUXILIARY_SET_INTERSECTOR_HPP_
