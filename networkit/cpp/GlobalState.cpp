
#include <networkit/GlobalState.hpp>

namespace NetworKit {

namespace GlobalState {

namespace {
std::atomic<uint64_t> seedValue{0};
std::atomic<uint64_t> globalSeedGeneration{0};
std::atomic<bool> seedUseThreadId{false};
} // namespace

void setSeed(uint64_t seed) {
    seedValue.store(seed, std::memory_order_relaxed);
}
uint64_t getSeed() {
    return seedValue.load(std::memory_order_relaxed);
}

uint64_t getGlobalSeed() {
    return globalSeedGeneration.load(std::memory_order_acquire);
}
void incGlobalSeed() {
    globalSeedGeneration.fetch_add(1, std::memory_order_release);
}

bool getSeedUseThreadId() {
    return seedUseThreadId.load(std::memory_order_relaxed);
}
void setSeedUseThreadId(bool useThreadId) {
    seedUseThreadId.store(useThreadId, std::memory_order_relaxed);
}

namespace {
bool printTime = false;
bool printLocation = false;

#ifdef NETWORKIT_QUIET_LOGGING
::Aux::Log::LogLevel loglevel = ::Aux::Log::LogLevel::quiet;
#else
::Aux::Log::LogLevel loglevel = ::Aux::Log::LogLevel::info;
#endif

std::ofstream logfile;

std::atomic_bool logfileIsOpen{false};
std::mutex logfileMutex;
} // namespace

Aux::Log::LogLevel getLogLevel() {
    return loglevel;
}

void setLogLevel(Aux::Log::LogLevel p) {
    loglevel = p;
}

bool getPrintTime() {
    return printTime;
}
void setPrintTime(bool b) {
    printTime = b;
}

bool getPrintLocation() {
    return printLocation;
}
void setPrintLocation(bool b) {
    printLocation = b;
}

std::ofstream &getLogFile() {
    return logfile;
}
void setLogfile(const std::string &filename) {
    std::lock_guard<std::mutex> guard{logfileMutex};
    if (logfile.is_open()) {
        logfile.close();
    }
    if (filename.empty()) {
        logfile.open(filename, std::ios_base::out | std::ios_base::app);
        logfileIsOpen = logfile.is_open();
    } else {
        logfileIsOpen = false;
    }
}

std::mutex &getLogFileMutex() {
    return logfileMutex;
}

bool getLogFileIsOpen() {
    return logfileIsOpen;
}

namespace {
std::atomic<bool> receivedSIGINT(false);
std::atomic<bool> rootSet(false);
std::atomic<Aux::SignalHandler *> root(0);
std::atomic<void (*)(int)> prevHandler;
} // namespace

bool getReceivedSIGINT() {
    return receivedSIGINT;
}
void setReceivedSIGINT(bool isReceivedSIGINT) {
    receivedSIGINT = isReceivedSIGINT;
}

bool getRootSet() {
    return rootSet;
}
void setRootSet(bool isRootSet) {
    rootSet = isRootSet;
}

auto getPrevHandler() -> void (*)(int) {
    return prevHandler;
}
void setPrevHandler(void (*handler)(int)) {
    prevHandler = handler;
}

Aux::SignalHandler *getRoot() {
    return root;
}
void setRoot(Aux::SignalHandler *rootHandler) {
    root = rootHandler;
}

} // namespace GlobalState

} // namespace NetworKit