// no-networkit-format
/*
 * Timer.h
 *
 *  Created on: 14.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_AUXILIARY_TIMER_HPP_
#define NETWORKIT_AUXILIARY_TIMER_HPP_

#include <chrono>
#include <cstdint>
#include <string>

#include <networkit/auxiliary/Log.hpp>

namespace Aux {

/**
 * A timer for running time measurements.
 */
class Timer {
public:
    #ifdef __MIC__
        using my_steady_clock = std::chrono::monotonic_clock;
    #else
        using my_steady_clock = std::chrono::steady_clock;
    #endif // __MIC__

    Timer() = default;

    /**
     * Start the clock.
     * Returns the time at which the instance was started.
     */
    my_steady_clock::time_point start() noexcept;

    /**
     * Stops the clock permanently for the instance of the Timer.
     * Returns the time at which the instance was stopped.
     */
    my_steady_clock::time_point stop() noexcept;

    /**
     * Returns a chrono::duration since the Timer was started.
     * If stop() was called, the duration is between the start() and stop()
     * calls is returned.
     */
    std::chrono::duration<uint64_t, std::milli> elapsed() const noexcept;

    /**
     * Returns the number of milliseconds since the Timer was started.
     * If stop() was called, the duration is between the start() and stop()
     * calls is returned.
     */
    uint64_t elapsedMilliseconds() const noexcept;

    /**
     * Returns the number of microseconds since the Timer was started.
     * If stop() was called, the duration is between the start() and stop()
     * calls is returned.
     */
    uint64_t elapsedMicroseconds() const noexcept;

    /**
     * Returns the number of nanoseconds since the Timer was started.
     * If stop() was called, the duration is between the start() and stop()
     * calls is returned.
     */
    uint64_t elapsedNanoseconds() const noexcept;

    /**
     * Returns the time at which the instance was started.
     */
    my_steady_clock::time_point startTime() const noexcept;

    /**
     * Returns the time at which the instance was stopped.
     */
    my_steady_clock::time_point stopTime() const noexcept;

    /**
     * Returns a human-readable representation including the elapsed time and unit.
     */
    std::string elapsedTag() const;

protected:
    bool running{false};                   //!< true if timer has been started and not stopped after that
    my_steady_clock::time_point started;   //!< time at which timer has been started
    my_steady_clock::time_point stopped;   //!< time at which timer has been stopped

    /// If running returns now, otherwise the stop time
    my_steady_clock::time_point stopTimeOrNow() const noexcept;
};

/**
 * A timer for running time measurements.
 * Same as Timer but automatically starts on construction.
 */
class StartedTimer : public Timer {
public:
    StartedTimer() : Timer() {
        start();
    }
};

/**
 * A timer for running time measurements within a scope.
 *
 * Same as Timer but automatically starts on construction and report on destruction.
 *
 * @code
 * {
 *    Aux::ScopedTimer("Algorithm A"); // WRONG; will directly report without measuring the scope
 *    Aux::ScopedTimer someName("Algorithm B"); // OK: the named instance is valid until the end of the scope
 *    // some expensive operations
 *
 *
 * } // <- Report time since creation of the timer object via Aux::Log
 * @endcode
 *
 * @warning
 * As the timer exploits RAII, it requires a named instance to extend its life time until the
 * end of the current scope. In the example above "someName" is used. If the time measured seems
 * too short, make you created a named instances.
 */
class LoggingTimer : public Timer {
public:
    /**
     * @param label The label printed out next to the measurement
     * @param level Only measure if logging at the provided LogLevel is enabled during runtime.
     *              If logging at the given level is disable the Timer is very cheap to construct.
     *              If logging is disabled during construction or destruction, no message is shown.
     */
    explicit LoggingTimer(const std::string& label = "", Aux::Log::LogLevel level = Aux::Log::LogLevel::debug);
    ~LoggingTimer();

private:
    Aux::Log::LogLevel level;
    std::string label;
};

} /* namespace Aux */

#endif // NETWORKIT_AUXILIARY_TIMER_HPP_

