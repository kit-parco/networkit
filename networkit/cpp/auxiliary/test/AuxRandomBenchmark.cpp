// no-networkit-format
#include <gtest/gtest.h>

#include <tlx/unused.hpp>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/Timer.hpp>

#include <atomic>
#include <random>
#include <omp.h>

namespace NetworKit {

class AuxRandomBenchmark: public testing::Test {
public:
    template<typename F>
    static double measure(F f, const size_t iterations = 50000000) {
        Aux::Timer timer;
        timer.start();

        for (size_t i = 0; i < iterations; i++) {
            f();
        }

        timer.stop();
        const auto ms = timer.elapsedMilliseconds();
        return (1.0e6 * ms / iterations);
    }

    template<typename F>
    static double measureParallel(F f) {
        // TODO: replace with google benchmark infrastructure
        std::atomic<uint64_t> atime{0}; // this is a very dirty hack, but atomic float-points are not fully support by standard
        std::atomic<int> num_threads;

        #pragma omp parallel
        {
            const double local_time = f();
            num_threads.store(omp_get_num_threads());
            atime.fetch_add(static_cast<uint64_t>(1e6 * local_time), std::memory_order_relaxed);
        }

        return 1e-6 * atime / num_threads;
    }

    template <typename T>
    static void DoNotOptimize(T&& x) {
        // TODO: replace with google benchmark infrastructure
        volatile auto dummy = x;
        tlx::unused(dummy);
    }
};

TEST_F(AuxRandomBenchmark, benchmarkInteger) {
    uint64_t tmp = 0;
    auto atime = measure([&] {
        tmp += Aux::Random::integer();
    });
    DoNotOptimize(tmp);
    std::cout << "Average time of operation: " << atime << "ns\n";
}

TEST_F(AuxRandomBenchmark, benchmarkLocalDistrInteger) {
    uint64_t tmp = 0;

    auto &prng = Aux::Random::getURNG();
    std::uniform_int_distribution<uint64_t> distr;

    auto atime = measure([&] {
        tmp += distr(prng);
    });

    DoNotOptimize(tmp);
    std::cout << "Average time of operation: " << atime << "ns\n";
}

TEST_F(AuxRandomBenchmark, benchmarkIntegerParallel) {
    auto atime = measureParallel([] {
        uint64_t tmp = 0;
        auto local_time = measure([&] {
            tmp += Aux::Random::integer();
        });
        DoNotOptimize(tmp);

        return local_time;
    });

    std::cout << "Average time of operation: " << atime << "ns\n";
}

TEST_F(AuxRandomBenchmark, benchmarkLocalDistrIntegerParallel) {
    auto atime = measureParallel([] {
        auto &prng = Aux::Random::getURNG();
        std::uniform_int_distribution<uint64_t> distr;

        uint64_t tmp = 0;
        auto local_time = measure([&] {
            tmp += distr(prng);
        });
        DoNotOptimize(tmp);

        return local_time;
    });

    std::cout << "Average time of operation: " << atime << "ns\n";
}


TEST_F(AuxRandomBenchmark, benchmarkProb) {
    double tmp = 0.0;
    auto atime = measure([&] {
        tmp += Aux::Random::probability();
    });
    DoNotOptimize(tmp);
    std::cout << "Average time of operation: " << atime << "ns\n";
}

TEST_F(AuxRandomBenchmark, benchmarkLocalDistrProb) {
    double tmp = 0.0;

    auto &prng = Aux::Random::getURNG();
    std::uniform_real_distribution<double> distr{0.0, 1.0};

    auto atime = measure([&] {
        tmp += distr(prng);
    });

    DoNotOptimize(tmp);
    std::cout << "Average time of operation: " << atime << "ns\n";
}

TEST_F(AuxRandomBenchmark, benchmarkProbParallel) {
    auto atime = measureParallel([] {
        double tmp = 0.0;
        auto local_time =  measure([&] {
            tmp += Aux::Random::probability();
        });
        DoNotOptimize(tmp);

        return local_time;
    });

    std::cout << "Average time of operation: " << atime << "ns\n";
}

TEST_F(AuxRandomBenchmark, benchmarkLocalDistrProbParallel) {
    auto atime = measureParallel([] {
        auto &prng = Aux::Random::getURNG();
        std::uniform_real_distribution<double> distr{0.0, 1.0};

        double tmp = 0.0;
        auto local_time = measure([&] {
            tmp += distr(prng);
        });
        DoNotOptimize(tmp);

        return local_time;
    });

    std::cout << "Average time of operation: " << atime << "ns\n";
}

} // ! namespace NetworKit
