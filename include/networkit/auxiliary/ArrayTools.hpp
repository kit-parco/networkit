#ifndef NETWORKIT_AUXILIARY_ARRAY_TOOLS_HPP_
#define NETWORKIT_AUXILIARY_ARRAY_TOOLS_HPP_

#include <cmath>
#include <utility>
#include <vector>

#include <networkit/Globals.hpp>

#include <tlx/math/integer_log2.hpp>

namespace Aux {
namespace ArrayTools {

template <class ArrayIt, class PermIt>
void applyPermutation(ArrayIt first, ArrayIt last, PermIt permFirst) {

    using PermType = typename std::iterator_traits<PermIt>::value_type;
    static_assert(std::is_integral<PermType>::value, "Elements of permutation must be integral.");

    const size_t arraySize =
        static_cast<typename std::iterator_traits<ArrayIt>::difference_type>(last - first);

    const size_t usedBitsInPerm = tlx::integer_log2_ceil(arraySize);
    // Avoid to use the sign bit if signed integral
    const size_t bitsInPerm = sizeof(PermType) * 8 - (std::is_signed<PermType>::value ? 1 : 0);

    if (bitsInPerm == usedBitsInPerm) {
        // All bits in permutation are needed, we mark swapped elements with a bool vector
        std::vector<bool> swapped(arraySize);
        for (size_t i = 0; i < arraySize; ++i) {
            if (swapped[i])
                continue;

            size_t cur = i;
            swapped[cur] = true;
            auto tmp = std::move(first[cur]);

            while (static_cast<size_t>(permFirst[cur]) != i) {
                first[cur] = std::move(first[permFirst[cur]]);
                cur = permFirst[cur];
                swapped[cur] = true;
            }

            first[cur] = std::move(tmp);
        }
    } else {
        // Not all bits in the permutation are needed, we use the last (or, if signed integral,
        // second to last) bit in the permutation to mark swapped elements
        const auto mask = PermType{1} << (bitsInPerm - (std::is_signed<PermType>::value ? 2 : 1));

        const auto swapped = [&permFirst, mask = mask](size_t i) -> bool {
            return (permFirst[i] & mask) != 0;
        };

        const auto markSwapped = [&permFirst, mask = mask](size_t i) -> void {
            permFirst[i] |= mask;
        };

        for (size_t i = 0; i < arraySize; ++i) {
            if (swapped(i))
                continue;

            size_t cur = i;
            markSwapped(cur);
            auto tmp = std::move(first[cur]);

            while (static_cast<size_t>(permFirst[cur] & ~mask) != i) {
                first[cur] = std::move(first[permFirst[cur] & ~mask]);
                cur = permFirst[cur] & ~mask;
                markSwapped(cur);
            }

            first[cur] = std::move(tmp);
        }

        for (size_t i = 0; i < arraySize; ++i)
            permFirst[i] &= ~mask;
    }
}

} // namespace ArrayTools
} // namespace Aux

#endif // NETWORKIT_AUXILIARY_ARRAY_TOOLS_HPP_
