/*
 * UniformRandomSelector.cpp
 *
 * Created: 2019-01-21
 * Author: Armin Wiebigke
 */

#include <random>

#include <networkit/auxiliary/UniformRandomSelector.hpp>
#include <networkit/auxiliary/Random.hpp>

namespace Aux {

UniformRandomSelector::UniformRandomSelector() : counter(0) {
}

bool UniformRandomSelector::addElement() {
    ++counter;
    std::uniform_int_distribution<NetworKit::count> dist(0, counter - 1);
    return (dist(Random::getURNG()) == 0);
}

void UniformRandomSelector::reset() {
    counter = 0;
}

}