
#include <exception>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/base/Algorithm.hpp>

namespace NetworKit {
std::string Algorithm::toString() const {
    throw std::runtime_error("TODO: implement in subclass and return string representation");
}

bool Algorithm::isParallel() const {
    throw std::runtime_error("TODO: Implement in subclass");
    return false;
}

} // namespace NetworKit
