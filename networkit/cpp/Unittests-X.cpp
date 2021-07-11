// no-networkit-format
/*
 * Unittests-X.cpp
 *
 *  Created on: 2012
 *      Author: Christian Staudt <christian.staudt@kit.edu>,
 *              Henning Meyerhenke <henning.meyerhenke@kit.edu>,
 *              Manuel Penschuck <networkit@manuel.jetzt>
 */

#include <algorithm>
#include <iostream>
#include <omp.h>
#include <utility>

#include <gtest/gtest.h>

#include <tlx/cmdline_parser.hpp>

#include <networkit/GlobalState.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallelism.hpp>

struct Options {
    std::string loglevel = "ERROR";
    bool sourceLocation{false};
    unsigned numThreads{0};

    bool modeTests{false};
    bool modeDebug{false};
    bool modeBenchmarks{false};
    bool modeRunnable{false};

    bool parse(int argc, char* argv[]) {
        tlx::CmdlineParser parser;

        parser.add_bool('t', "tests",      modeTests,      "Run unit tests");
        parser.add_bool('r', "run",        modeRunnable,   "Run unit tests which don't use assertions");
        parser.add_bool('d', "debug",      modeDebug,      "Run tests to debug some algorithms");
        parser.add_bool('b', "benchmarks", modeBenchmarks, "Run benchmarks");

        parser.add_unsigned("threads",     numThreads,     "set the maximum number of threads; 0 (=default) uses OMP default");
        parser.add_string("loglevel",      loglevel,       "set the log level (TRACE|DEBUG|INFO|WARN|ERROR|FATAL)");
        parser.add_bool("srcloc",          sourceLocation, "print source location of log messages");

        if (!parser.process(argc, argv, std::cerr))
            return false;

        if (modeTests + modeDebug + modeBenchmarks + modeRunnable > 1) {
            std::cerr << "Select at most one of -t, -r, -d, -b\n";
            return false;
        }

        return true;
    }
};

int main(int argc, char *argv[]) {
    std::cout << "*** NetworKit Unit Tests ***\n";

    ::testing::InitGoogleTest(&argc, argv);
    Options options;
    if (!options.parse(argc, argv)) {
        return -1;
    }

    // Configure logging
    Aux::Log::setLogLevel(options.loglevel);
    NetworKit::GlobalState::setPrintLocation(options.sourceLocation);
    std::cout << "Loglevel: " << Aux::Log::getLogLevel() << "\n";
#ifdef NETWORKIT_RELEASE_LOGGING
    if (options.loglevel == "TRACE" || options.loglevel == "DEBUG")
        std::cout << "WARNING: TRACE and DEBUG messages are missing"
                     " in NETWORKIT_RELEASE_LOGGING builds" << std::endl;
#endif

    // Configure parallelism
    {
        if (options.numThreads) {
            Aux::setNumberOfThreads(options.numThreads);
        }
        std::cout << "Max. number of threads: " << Aux::getMaxNumberOfThreads() << "\n";
    }

    // Configure test filter
    {
        auto setFilter = [&](const char *filter) {
            ::testing::GTEST_FLAG(filter) = filter;
        };

        if (options.modeTests) {
            setFilter("*Test.test*");
        } else if (options.modeDebug) {
            setFilter("*Test.debug*");
        } else if (options.modeBenchmarks) {
            setFilter("*Benchmark*");
        } else if (options.modeRunnable) {
            setFilter("*Test.run*");
        }
    }

    INFO("=== starting unit tests ===");

    return RUN_ALL_TESTS();
}
