
#ifndef NETWORKIT_GLOBAL_STATE_HPP_
#define NETWORKIT_GLOBAL_STATE_HPP_

#if defined(NETWORKIT_BUILDING_STATELIB)
#if defined(WIN32)
#define NETWORKIT_EXPORT __declspec(dllexport)
#else
#define NETWORKIT_EXPORT
#endif // WIN32
#else
#if defined(WIN32)
#define NETWORKIT_EXPORT __declspec(dllimport)
#else
#define NETWORKIT_EXPORT
#endif // WIN32
#endif // NETWORKIT_BUILDING_STATELIB

#include <atomic>
#include <csignal>
#include <cstdint>
#include <fstream>
#include <mutex>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>

namespace NetworKit {

namespace GlobalState {

// Global states used for auxiliary/Random

// The global seed generation functions as an epoch counter
// for each time, the seed was updated. At the same time it
// releases the changes made to the seed values to other
// threads which always access it by an aquire.
// seedValues and seedUseThreadId can hence be accesses
// relaxed.
// As long as globalSeedGeneration is zero, the getURNG()
// ignores these two variables.

NETWORKIT_EXPORT uint64_t getSeed();
NETWORKIT_EXPORT void setSeed(uint64_t seed);

NETWORKIT_EXPORT uint64_t getGlobalSeed();
NETWORKIT_EXPORT void incGlobalSeed();

NETWORKIT_EXPORT void setSeedUseThreadId(bool useThreadId);
NETWORKIT_EXPORT bool getSeedUseThreadId();

// Global states used for auxiliary/Log
NETWORKIT_EXPORT Aux::Log::LogLevel getLogLevel();
NETWORKIT_EXPORT void setLogLevel(Aux::Log::LogLevel p = Aux::Log::LogLevel::info);

NETWORKIT_EXPORT bool getPrintTime();
NETWORKIT_EXPORT void setPrintTime(bool b);

NETWORKIT_EXPORT bool getPrintLocation();
NETWORKIT_EXPORT void setPrintLocation(bool b);

NETWORKIT_EXPORT std::ofstream &getLogFile();
NETWORKIT_EXPORT void setLogfile(const std::string &filename);

NETWORKIT_EXPORT std::mutex &getLogFileMutex();

NETWORKIT_EXPORT bool getLogFileIsOpen();

// Global states used for auxiliary/SignalHandling
NETWORKIT_EXPORT bool getReceivedSIGINT();
NETWORKIT_EXPORT void setReceivedSIGINT(bool isReceivedSIGINT);

NETWORKIT_EXPORT bool getRootSet();
NETWORKIT_EXPORT void setRootSet(bool isRootSet);

NETWORKIT_EXPORT auto getPrevHandler() -> void (*)(int);
NETWORKIT_EXPORT void setPrevHandler(void (*handler)(int));

NETWORKIT_EXPORT Aux::SignalHandler *getRoot();
NETWORKIT_EXPORT void setRoot(Aux::SignalHandler *rootHandler);

} // namespace GlobalState

} // namespace NetworKit

#endif // NETWORKIT_GLOBAL_STATE_HPP_
