// no-networkit-format
/*
 * Random.cpp
 *
 *  Created on: 02.01.2014
 *      Author: FJW
 */

#include <cmath>
#include <omp.h>
#include <limits>
#include <atomic>

#include <networkit/GlobalState.hpp>
#include <networkit/auxiliary/Random.hpp>

namespace Aux {
namespace Random {

void setSeed(uint64_t seed, bool useThreadId) {
    NetworKit::GlobalState::setSeed(seed);
    NetworKit::GlobalState::setSeedUseThreadId(useThreadId);
    NetworKit::GlobalState::incGlobalSeed();
    getURNG(); // update local seed value
}

static uint64_t getSeed_(uint64_t globSeedGen) {
    if (!globSeedGen) {
        thread_local static std::random_device urng{};
        std::uniform_int_distribution<uint64_t> dist{};
        return dist(urng);

    } else if (NetworKit::GlobalState::getSeedUseThreadId()) {
        return NetworKit::GlobalState::getSeed() + omp_get_thread_num();

    } else {
        return NetworKit::GlobalState::getSeed();
    }
}

uint64_t getSeed() {
    return getSeed_(NetworKit::GlobalState::getGlobalSeed());
}

bool getUseThreadId() {
    return NetworKit::GlobalState::getSeedUseThreadId();
}

std::mt19937_64& getURNG() {
    thread_local static std::mt19937_64 generator{getSeed()};
    thread_local static uint64_t localSeedGeneration{0};

    auto globSeedGen = NetworKit::GlobalState::getGlobalSeed();

    if (localSeedGeneration != globSeedGen) {
        generator.seed(getSeed_(globSeedGen));
        localSeedGeneration = globSeedGen;
    }
    return generator;
}

uint64_t integer() {
    thread_local static std::uniform_int_distribution<uint64_t> dist{};
    return dist(getURNG());
}
uint64_t integer(uint64_t upperBound) {
    std::uniform_int_distribution<uint64_t> dist{0, upperBound};
    return dist(getURNG());
}
uint64_t integer(uint64_t lowerBound, uint64_t upperBound) {
    std::uniform_int_distribution<uint64_t> dist{lowerBound, upperBound};
    return dist(getURNG());
}

double real() {
    thread_local static std::uniform_real_distribution<double> dist{};
    return dist(getURNG());
}
double real(double upperBound) {
    std::uniform_real_distribution<double> dist{0.0, upperBound};
    return dist(getURNG());
}
double real(double lowerBound, double upperBound) {
    std::uniform_real_distribution<double> dist{lowerBound, upperBound};
    return dist(getURNG());
}

double probability() {
    thread_local static std::uniform_real_distribution<double> dist{0.0, std::nexttoward(1.0, 2.0)};
    return dist(getURNG());
}

std::size_t index(std::size_t max) {
    assert(max > 0 && "There have to be valid indexes");
    std::uniform_int_distribution<std::size_t> dist{0, max - 1};
    return dist(getURNG());
}

} // namespace Random
} // namespace Aux
