/*
 * UniformRandomSelector.hpp
 *
 * Created: 2019-01-21
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_AUXILIARY_UNIFORM_RANDOM_SELECTOR_HPP_
#define NETWORKIT_AUXILIARY_UNIFORM_RANDOM_SELECTOR_HPP_

#include <random>

#include <networkit/Globals.hpp>

namespace Aux {

class UniformRandomSelector {
public:
    UniformRandomSelector();

    bool addElement();

    void reset();

private:
    NetworKit::count counter;
};

}

#endif // NETWORKIT_AUXILIARY_UNIFORM_RANDOM_SELECTOR_HPP_
