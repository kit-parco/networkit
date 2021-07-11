#ifndef NETWORKIT_AUXILIARY_HASH_UTILS_HPP_
#define NETWORKIT_AUXILIARY_HASH_UTILS_HPP_

#include <utility>

namespace Aux {

template <typename T>
void hashCombine(std::size_t &seed, T const &v) {
    // The Boost Software License, Version 1.0 applies to this function.
    // See https://www.boost.org/LICENSE_1_0.txt
    // https://www.boost.org/doc/libs/1_75_0/doc/html/hash/reference.html#boost.hash_combine
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct PairHash {
    template <typename A, typename B>
    std::size_t operator()(const std::pair<A, B> &pair) const {
        std::size_t seed = 0;
        hashCombine(seed, pair.first);
        hashCombine(seed, pair.second);
        return seed;
    }
};

} // namespace Aux

#endif // NETWORKIT_AUXILIARY_HASH_UTILS_HPP_
