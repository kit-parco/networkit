// no-networkit-format
/*
 * Header which provides parallel STL implementations when available.
 */
#ifndef NETWORKIT_AUXILIARY_PARALLEL_HPP_
#define NETWORKIT_AUXILIARY_PARALLEL_HPP_

#include <atomic>
#include <functional>

#if ((defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))) && defined _OPENMP
#include <parallel/algorithm>
#else
#define NETWORKIT_NO_PARALLEL_STL
#include <algorithm>
#endif

namespace Aux {
    namespace Parallel {
        #ifdef NETWORKIT_NO_PARALLEL_STL
        using std::sort;
        using std::max_element;
        #else
        using __gnu_parallel::sort;
        using __gnu_parallel::max_element;
        #endif // NETWORKIT_NO_PARALLEL_STL

        template <typename ValueType, typename Comp>
        void atomic_set(std::atomic<ValueType> &target, const ValueType &input, Comp shallSet) {
            // load current value
            auto curValue = target.load(std::memory_order_relaxed);

            do {
                if (!shallSet(input, curValue)) break; // skip if shall not be set (anymore)
            // set new value unless current value has been changed in the meantime, if current value has changed load it (so we can compare again)
            } while (!target.compare_exchange_weak(curValue, input, std::memory_order_release, std::memory_order_relaxed));
        }

        template <typename ValueType>
        void atomic_max(std::atomic<ValueType> &target, const ValueType &input) {
            atomic_set(target, input, std::greater<ValueType>());
        }

        template <typename ValueType>
        void atomic_min(std::atomic<ValueType> &target, const ValueType &input) {
            atomic_set(target, input, std::less<ValueType>());
        }

    }
}

#endif // NETWORKIT_AUXILIARY_PARALLEL_HPP_
