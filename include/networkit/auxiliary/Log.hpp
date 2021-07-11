// no-networkit-format
#ifndef NETWORKIT_AUXILIARY_LOG_HPP_
#define NETWORKIT_AUXILIARY_LOG_HPP_

#include <sstream>
#include <string>

#include <networkit/auxiliary/StringBuilder.hpp>

#include <tlx/define/deprecated.hpp>

#ifdef _MSC_VER
    #define NETWORKT_PRETTY_FUNCTION __FUNCSIG__
#else
    #define NETWORKT_PRETTY_FUNCTION __PRETTY_FUNCTION__
#endif // _MSC_VER

/// Logging without format string
#define LOG_AT(level, ...) ::Aux::Log::log({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__}, level, __VA_ARGS__)
#define FATAL(...) LOG_AT(::Aux::Log::LogLevel::fatal, __VA_ARGS__)
#define ERROR(...) LOG_AT(::Aux::Log::LogLevel::error, __VA_ARGS__)
#define WARN(...)  LOG_AT(::Aux::Log::LogLevel::warn,  __VA_ARGS__)
#define INFO(...)  LOG_AT(::Aux::Log::LogLevel::info,  __VA_ARGS__)

/// Logging with format string
#define LOG_ATF(level, ...) ::Aux::Log::logF({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},level, __VA_ARGS__)
#define FATALF(...) LOG_ATF(::Aux::Log::LogLevel::fatal, __VA_ARGS__)
#define ERRORF(...) LOG_ATF(::Aux::Log::LogLevel::error, __VA_ARGS__)
#define WARNF(...)  LOG_ATF(::Aux::Log::LogLevel::warn,  __VA_ARGS__)
#define INFOF(...)  LOG_ATF(::Aux::Log::LogLevel::info,  __VA_ARGS__)

// DEBUG and TRACE are no-ops if NETWORKIT_RELEASE_LOGGING is defined.
#if defined(NETWORKIT_RELEASE_LOGGING)
#   define NETWORKIT_LOG_DUMMY(...) do {} while(false)
#	define DEBUG(...)  NETWORKIT_LOG_DUMMY(__VA_ARGS__)
#	define DEBUGF(...) NETWORKIT_LOG_DUMMY(__VA_ARGS__)
#	define TRACE(...)  NETWORKIT_LOG_DUMMY(__VA_ARGS__)
#	define TRACEF(...) NETWORKIT_LOG_DUMMY(__VA_ARGS__)
#	define TRACEPOINT  NETWORKIT_LOG_DUMMY()
#else
#	define DEBUG(...) LOG_AT(::Aux::Log::LogLevel::debug, __VA_ARGS__)
#	define TRACE(...) LOG_AT(::Aux::Log::LogLevel::trace, __VA_ARGS__)

#	define DEBUGF(...) LOG_ATF(::Aux::Log::LogLevel::debug, __VA_ARGS__)
#	define TRACEF(...) LOG_ATF(::Aux::Log::LogLevel::trace, __VA_ARGS__)

#	define TRACEPOINT LOG_AT(::Aux::Log::LogLevel::trace, "tracepoint")
#endif // defined(NETWORKIT_RELEASE_LOGGING)

namespace Aux { namespace Log {

struct Location {
    const char* file;
    const char* function;
    const int line;
};

enum class LogLevel {
    trace,
    debug,
    info,
    warn,
    error,
    fatal,
    quiet // Emits no log messages at all.
};

/**
 * Accept loglevel as string and set.
 * @param logLevel as string
 */
void setLogLevel(const std::string &logLevel);

/**
 * @return current loglevel as string
 */
std::string getLogLevel();

namespace Settings {

LogLevel TLX_DEPRECATED(getLogLevel());

void TLX_DEPRECATED(setLogLevel(LogLevel p = LogLevel::info));

bool TLX_DEPRECATED(getPrintTime());

void TLX_DEPRECATED(setPrintTime(bool b));

bool TLX_DEPRECATED(getPrintLocation());

void TLX_DEPRECATED(setPrintLocation(bool b));

void TLX_DEPRECATED(setLogfile(const std::string &filename));
} // namespace Settings

namespace Impl {
void log(const Location &loc, LogLevel p, const std::string &msg);
} //namespace impl

///! Returns true iff logging at the provided level is currently activated
bool isLogLevelEnabled(LogLevel p) noexcept;

template<typename...T>
void log(const Location &loc, LogLevel p, const T &...args) {
    if(!isLogLevelEnabled(p))
        return;

    std::stringstream stream;
    printToStream(stream, args...);
    Impl::log(loc, p, stream.str());
}

template<typename...T>
void logF(const Location &loc, LogLevel p, const std::string &format, const T &...args) {
    if(!isLogLevelEnabled(p))
        return;

    std::stringstream stream;
    printToStreamF(stream, format, args...);
    Impl::log(loc, p, stream.str());
}

}} // namespace Aux::Log

#endif // NETWORKIT_AUXILIARY_LOG_HPP_
