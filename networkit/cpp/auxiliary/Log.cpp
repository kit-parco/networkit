// no-networkit-format
#include <atomic>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ios>
#include <mutex>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/GlobalState.hpp>

namespace Aux { namespace Log {

void setLogLevel(const std::string &logLevel) {
    if (logLevel == "TRACE") {
        NetworKit::GlobalState::setLogLevel(LogLevel::trace);
    } else if (logLevel == "DEBUG") {
        NetworKit::GlobalState::setLogLevel(LogLevel::debug);
    } else if (logLevel == "INFO") {
        NetworKit::GlobalState::setLogLevel(LogLevel::info);
    } else if (logLevel == "WARN") {
        NetworKit::GlobalState::setLogLevel(LogLevel::warn);
    } else if (logLevel == "ERROR") {
        NetworKit::GlobalState::setLogLevel(LogLevel::error);		
    } else if (logLevel == "FATAL") {
        NetworKit::GlobalState::setLogLevel(LogLevel::fatal);
    } else if (logLevel == "QUIET") {
        NetworKit::GlobalState::setLogLevel(LogLevel::quiet);
    } else {
        throw std::runtime_error("unknown loglevel");
    }
}

std::string getLogLevel() {
    LogLevel current = NetworKit::GlobalState::getLogLevel();
    switch (current) {
    case LogLevel::trace:
        return "TRACE";
    case LogLevel::debug:
        return "DEBUG";
    case LogLevel::info:
        return "INFO";
    case LogLevel::warn:
        return "WARN";
    case LogLevel::error:
        return "ERROR";
    case LogLevel::fatal:
        return "FATAL";
    case LogLevel::quiet:
        return "QUIET";
    default:
        throw std::logic_error{"invalid loglevel in getLogLevel()"};
    }
}

namespace Settings {

LogLevel getLogLevel() {
    return NetworKit::GlobalState::getLogLevel();
}

void setLogLevel(LogLevel p) {
    NetworKit::GlobalState::setLogLevel(p);    
}

bool getPrintTime() {
    return NetworKit::GlobalState::getPrintTime();  
}

void setPrintTime(bool b) {
    NetworKit::GlobalState::setPrintTime(b);    
}

bool getPrintLocation() {
    return NetworKit::GlobalState::getPrintLocation();
}

void setPrintLocation(bool b) {
    NetworKit::GlobalState::setPrintLocation(b);   
}

void setLogfile(const std::string &filename) {
    NetworKit::GlobalState::setLogfile(filename);    
} 
}

void printLogLevel(std::ostream &stream, LogLevel p) {
    switch(p) {
        case LogLevel::fatal:
            stream << "[FATAL]"; break;
        case LogLevel::error:
            stream << "[ERROR]"; break;
        case LogLevel::warn:
            stream << "[WARN ]"; break;
        case LogLevel::info:
            stream << "[INFO ]"; break;
        case LogLevel::debug:
            stream << "[DEBUG]"; break;
        case LogLevel::trace:
            stream << "[TRACE]"; break;
        default:
            break;
    }
}

bool isLogLevelEnabled(LogLevel p) noexcept {
    return p >= NetworKit::GlobalState::getLogLevel();
}

void printTime(std::ostream &stream,
        const std::chrono::time_point<std::chrono::system_clock> &timePoint) {
    stream << "[" << timePoint.time_since_epoch().count() << "]";
}

void printLocation(std::ostream &stream, const Location &loc) {
    stream << "[" << loc.file << ", " << loc.line << ": " << loc.function << "]";
}

std::tuple<std::string, std::string> getTerminalFormat(LogLevel p) {
    switch (p) {
        case LogLevel::fatal:
            return std::make_tuple("\033[1;31m", "\033[0m");
        case LogLevel::error:
            return std::make_tuple("\033[31m", "\033[0m");
        case LogLevel::warn:
            return std::make_tuple("\033[33m", "\033[0m");
        case LogLevel::info:
        case LogLevel::debug:
        case LogLevel::trace:
            return std::make_tuple("", "");
        default:
            // this only exists to silence a warning:
            // TODO: consider replacing it with __builtin_unreachable();
            throw std::logic_error{"invalid loglevel. This should NEVER happen"};
    }
}

static void logToTerminal(const Location &loc, LogLevel p,
        const std::chrono::time_point<std::chrono::system_clock> &timePoint,
        const std::string &msg) {
    std::stringstream stream;
    
    if(NetworKit::GlobalState::getPrintTime()) {
        printTime(stream, timePoint);
    }
    
    std::string termFormatOpen, termFormatClose;
    std::tie(termFormatOpen, termFormatClose) = getTerminalFormat(p);
    
    stream << termFormatOpen;
    printLogLevel(stream, p);
    stream <<termFormatClose;
    
    if(NetworKit::GlobalState::getPrintLocation()) {
        printLocation(stream, loc);
    }
    
    stream << ": ";
    
    stream << termFormatOpen;
    stream << msg;
    stream << termFormatClose;
    
    stream.put('\n');
    
    static std::mutex cerr_mutex;
    {
        std::lock_guard<std::mutex> guard{cerr_mutex};
        std::cerr << stream.str();
    }
}

static void logToFile(const Location &loc, LogLevel p,
        const std::chrono::time_point<std::chrono::system_clock> &timePoint,
        const std::string &msg) {
    if(!NetworKit::GlobalState::getLogFileIsOpen()) {
        return;
    }
    std::stringstream stream;
    printTime(stream, timePoint);
    stream << ' ';
    printLogLevel(stream, p);
    
    if(NetworKit::GlobalState::getPrintLocation()) {
        stream << ' ';
        printLocation(stream, loc);
    }
    
    stream << ": " << msg << '\n';
    {
        std::lock_guard<std::mutex> guard{NetworKit::GlobalState::getLogFileMutex()};
        if(!NetworKit::GlobalState::getLogFileIsOpen()) {
            return;
        }
        NetworKit::GlobalState::getLogFile() << stream.str() << std::flush;
    }
}

namespace Impl {

void log(const Location &loc, LogLevel p, const std::string &msg) {
    auto time =std::chrono::system_clock::now();
    
    logToTerminal(loc, p, time, msg);
    logToFile(loc, p, time, msg);
}

} // namespace impl

}} // namespace Aux::Log
