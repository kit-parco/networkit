// no-networkit-format
/*
 * Timer.cpp
 *
 *  Created on: 14.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <iostream>
#include <sstream>

#include <networkit/auxiliary/Timer.hpp>

namespace Aux {

Timer::my_steady_clock::time_point Timer::start() noexcept {
    started = my_steady_clock::now();
    running = true;
    return started;
}

Timer::my_steady_clock::time_point Timer::stop() noexcept {
    stopped = my_steady_clock::now();
    running = false;
    return stopped;
}

std::chrono::duration<uint64_t, std::milli> Timer::elapsed() const noexcept {
    return std::chrono::duration_cast<std::chrono::duration<uint64_t, std::milli> >(stopTimeOrNow() - started);
}

Timer::my_steady_clock::time_point Timer::startTime() const noexcept {
    return started;
}

Timer::my_steady_clock::time_point Timer::stopTime() const noexcept {
    return stopped;
}

uint64_t Timer::elapsedMilliseconds() const noexcept {
    return elapsed().count();
}

uint64_t Timer::elapsedMicroseconds() const noexcept {
    return std::chrono::duration_cast<std::chrono::duration<uint64_t, std::micro>>(stopTimeOrNow() - started).count();
}

uint64_t Timer::elapsedNanoseconds() const noexcept {
    return std::chrono::duration_cast<std::chrono::duration<uint64_t, std::nano>>(stopTimeOrNow() - started).count();
}

std::string Timer::elapsedTag() const {
    std::stringstream s;
    s << "(" << elapsedMilliseconds() << " ms) ";
    return s.str();
}

Timer::my_steady_clock::time_point Timer::stopTimeOrNow() const noexcept {
    return running ? std::chrono::steady_clock::now() : stopped;
}

LoggingTimer::LoggingTimer(const std::string &label, Aux::Log::LogLevel level)
    : level(level)
{
    if (!Aux::Log::isLogLevelEnabled(level))
        return;

    this->label = label;
    start();
}

LoggingTimer::~LoggingTimer() {
    if (!running)
        return;

    std::stringstream ss;
    ss << "Timer ";

    if (!label.empty())
        ss << '"' << label << "\" ";

    ss << "ran for " << (static_cast<double>(elapsedMicroseconds()) * 1e-3) << " ms";

    LOG_AT(level, ss.str());
}

} /* namespace Aux */
