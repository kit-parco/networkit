/*
 * ParallelTimings.cpp
 *
 * Created on: 2019-11-05
 * Author: Armin Wiebigke
  */

#include <map>

#include <networkit/auxiliary/ParallelTimings.hpp>

namespace NetworKit {

void ParallelTimings::addTime(Aux::Timer &timer, const std::string &name) const {
    timer.stop();
    double elapsed = timer.elapsedNanoseconds();
    auto tid = omp_get_thread_num();
    timings[tid][name] += elapsed;
    timer.start();
}

ParallelTimings::ParallelTimings() : timings(omp_get_max_threads()) {

}

void ParallelTimings::addTimings(const std::unordered_map<std::string, double> &ts, const std::string &prefix) const {
    auto tid = omp_get_thread_num();
    for (auto &t : ts) {
        timings[tid][prefix + t.first] += t.second;
    }
}

std::unordered_map<std::string, double> ParallelTimings::getTimings() const {
    std::unordered_map<std::string, double> timingsSum;
    for (const auto& threadTimings : timings) {
        for (const auto& it : threadTimings) {
            timingsSum[it.first] += it.second;
        }
    }
    return timingsSum;
}

bool ParallelTimings::timingsEmpty() const {
    for (const auto& threadTimings : timings) {
        if (!threadTimings.empty())
            return false;
    }
    return true;
}

std::string ParallelTimings::timingsAsString() const {
    std::map<std::string, double> timingsSum;
    for (const auto& threadTimings : timings) {
        for (const auto& it : threadTimings) {
            timingsSum[it.first] += it.second;
        }
    }
    std::stringstream str;
    for (const auto& t : timingsSum) {
        str << t.first + ": " << t.second / 1e6 << "ms" << "\n";
    }
    return str.str();
}

} // namespace NetworKit

