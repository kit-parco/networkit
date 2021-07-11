// no-networkit-format
/*
 * Random.h
 *
 *  Created on: 02.01.2014
 *      Author: FJW
 */

#ifndef NETWORKIT_AUXILIARY_RANDOM_HPP_
#define NETWORKIT_AUXILIARY_RANDOM_HPP_

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>

#include <tlx/define/deprecated.hpp>

namespace Aux {

/**
 * Provides several functions for random-numbers.
 *
 * All functions are guaranteed to be thread-safe.
 */
namespace Random {

/**
 * Sets the random seed that is used globally.
 * @param seed The seed value
 * @param useThreadId If the thread id shall be added to the seed value
 */
void setSeed(uint64_t seed, bool useThreadId);

/**
 * @returns a high-quality random seed for an URNG.
 */
uint64_t getSeed();

/*
 * @returns whether or not the thread id is added to the random seed. 
 */
bool TLX_DEPRECATED(getUseThreadId());

/**
 * @returns a reference to a seeded URNG that is thread_local.
 */
std::mt19937_64& getURNG();


/**
 * @returns an integer distributed uniformly in an inclusive range;
 * @param upperBound the upper bound, default = UNINT64_T_MAX
 * @param lowerBound the lower bound, default = 0
 *
 * @warning Compared to obtaining a reference to a generator using
 * @ref getURNG() and then using a local std::uniform_int_distribution,
 * this method incurs a slow-down of up to 30%. Consider avoiding it
 * in hot sections.
 */
uint64_t integer();
uint64_t integer(uint64_t upperBound);
uint64_t integer(uint64_t lowerBound, uint64_t upperBound);

/**
 * @returns a double distributed uniformly in a half-open range: [lowerBound, upperBound)
 * @param upperBound default = 1.0
 * @param lowerBound default = 0.0
 *
 * @warning Compared to obtaining a reference to a generator using
 * @ref getURNG() and then using a local std::uniform_real_distribution,
 * this method incurs a slow-down of up to 30%. Consider avoiding it
 * in hot sections.
 */
double real();
double real(double upperBound);
double real(double lowerBound, double upperBound);

/**
 * @returns a double distributed uniformly in the range [0, 1]
 * @note this differs from real() in returning a value in a closed instead of a half-open range.
 *
 * @warning Compared to obtaining a reference to a generator using
 * @ref getURNG() and then using a local std::uniform_int_distribution,
 * this method incurs a slow-down of up to 30%. Consider avoiding it
 * in hot sections.
 */
double probability();

/**
 * @returns a size_t in the range [0, max - 1] that can for example be used to
 * access a random element of a container.
 */
std::size_t index(std::size_t max);

/**
 * @returns a uniform random choice from an indexable container of elements.
 */
template<typename Container>
typename Container::const_reference choice(const Container& container) {
    return container[index(container.size())];
}


/**
 * @returns a weighted random choice from a vector of elements with given weights.
 */
template <typename Element>
const Element& weightedChoice(const std::vector<std::pair<Element, double>>& weightedElements) {
    if (weightedElements.empty())
        throw std::runtime_error("Random::weightedChoice: input size equal to 0");
    double total = 0.0;
    for (const auto& entry : weightedElements) {
        assert(entry.second >= 0.0 && "This algorithm only works with non-negative weights");
        total += entry.second;
    }
    double r = Aux::Random::real(total);
    for (const auto& entry : weightedElements) {
        if (r < entry.second) {
            return entry.first;
        }
        r -= entry.second;
    }
    throw std::runtime_error("Random::weightedChoice: should never get here"); // should never get here
}

} // namespace Random
} // namespace Aux


#endif // NETWORKIT_AUXILIARY_RANDOM_HPP_
