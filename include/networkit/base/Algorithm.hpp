
#ifndef NETWORKIT_BASE_ALGORITHM_HPP_
#define NETWORKIT_BASE_ALGORITHM_HPP_

#include <stdexcept>
#include <string>

#include <tlx/define/deprecated.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

class Algorithm {
protected:
    bool hasRun = false;

public:
    Algorithm() = default;

    virtual ~Algorithm() = default;

    /**
     * The generic run method which calls runImpl() and takes care of setting @ref hasRun to the
     * appropriate value.
     */
    virtual void run() = 0;

    /**
     * Indicates whether an algorithm has completed computation or not.
     *
     * @return The value of @ref hasRun.
     */
    bool hasFinished() const noexcept { return hasRun; }

    /**
     * Assure that the algorithm has been run, throws a std::runtime_error otherwise.
     */
    void assureFinished() const {
        if (!hasRun)
            throw std::runtime_error("Error, run must be called first");
    }

    /**
     * Returns a string with the algorithm's name and its parameters, if there
     * are any. Subclasses should override it.
     *
     * @return The string representation of the algorithm.
     */
    virtual std::string TLX_DEPRECATED(toString() const);

    /**
     * @return True if algorithm can run multi-threaded.
     */
    virtual TLX_DEPRECATED(bool isParallel() const);
};

} // namespace NetworKit

#endif // NETWORKIT_BASE_ALGORITHM_HPP_
