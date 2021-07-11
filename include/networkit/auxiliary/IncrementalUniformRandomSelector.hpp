#ifndef NETWORKIT_AUXILIARY_INCREMENTAL_UNIFORM_RANDOM_SELECTOR_HPP_
#define NETWORKIT_AUXILIARY_INCREMENTAL_UNIFORM_RANDOM_SELECTOR_HPP_

#include <cstddef>

#include <networkit/auxiliary/Random.hpp>

namespace Aux {

/**
 * Select uniformly at random from a set of elements
 * that is discovered incrementally. Every time you
 * discover a new element, you tell the selector and
 * it will tell if you should keep the new element or
 * your old element. As more and more elements are
 * discovered, it will be less and less likely that
 * the new element is chosen in order to ensure that
 * the selection is uniformly at random. The process
 * can be reset when a new class of elements is found
 * that shall be used instead.
 */
class IncrementalUniformRandomSelector {
public:
    /**
     * Initialize the random select for one element.
     */
    IncrementalUniformRandomSelector() : counter(1){};

    /**
     * Add the next element.
     *
     * Note that this must not be called for the first element.
     *
     * @return If the new element shall be accepted/kept.
     */
    bool addElement() {
        ++counter;
        return (Random::real() < 1.0 / static_cast<double>(counter));
    }

    /**
     * Reset the process to one element again.
     */
    void reset() { counter = 1; }

private:
    size_t counter;
};

} // namespace Aux

#endif // NETWORKIT_AUXILIARY_INCREMENTAL_UNIFORM_RANDOM_SELECTOR_HPP_
