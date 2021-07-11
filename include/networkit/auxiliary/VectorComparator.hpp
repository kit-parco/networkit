// networkit-format

#ifndef NETWORKIT_AUXILIARY_VECTOR_COMPARATOR_HPP_
#define NETWORKIT_AUXILIARY_VECTOR_COMPARATOR_HPP_

#include <cstdint>
#include <vector>

namespace Aux {

// Implementation of vector-based comparators

template <class Type>
struct LessInVector {
    LessInVector(const std::vector<Type> &vec) : vec(&vec) {}
    bool operator()(uint64_t x, uint64_t y) const noexcept { return (*vec)[x] < (*vec)[y]; }

private:
    const std::vector<Type> *vec;
};

template <class Type>
struct GreaterInVector {
    GreaterInVector(const std::vector<Type> &vec) : vec(&vec) {}
    bool operator()(uint64_t x, uint64_t y) const noexcept { return (*vec)[x] > (*vec)[y]; }

private:
    const std::vector<Type> *vec;
};

} // namespace Aux

#endif // NETWORKIT_AUXILIARY_VECTOR_COMPARATOR_HPP_
