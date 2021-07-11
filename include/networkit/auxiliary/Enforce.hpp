// no-networkit-format
#ifndef NETWORKIT_AUXILIARY_ENFORCE_HPP_
#define NETWORKIT_AUXILIARY_ENFORCE_HPP_

#include <cassert>
#include <stdexcept>
#include <string>

namespace Aux {

/**
 * Enforces that @a b is true and throws an Exception otherwise.
 *
 * If provided, msg must not be null, otherwise the behavior of this
 * function is undefined.
 *
 * @param b Boolean value whose truthiness should be enforced
 * @param msg Message of the exception
 */
template<typename Exception = std::runtime_error>
inline void enforce(bool b, const char* msg = "") {
    assert(msg && "Message to enforce must not be nullptr");
    if (!b) {
        throw Exception{msg};
    }
}

/**
 * Enforces that @a b is true and throws an Exception otherwise.
 *
 * If provided, msg must not be null, otherwise the behavior of this
 * function is undefined.
 *
 * @param b Boolean value whose truthiness should be enforced
 * @param msg Message of the exception
 */
template<typename Exception = std::runtime_error>
inline void enforce(bool b, const std::string& msg) {
    enforce<Exception>(b, msg.c_str());
}

/**
 * Checks that the provided fstream is opened and throws an exception otherwise.
 *
 * @param stream File stream whose openness should be enforced
 */
template<typename Stream>
inline void enforceOpened(const Stream& stream) {
    enforce(stream.is_open());
}

/**
 * This namespace provides some Types with a static member-function `void enforce(bool)`
 * that may check wether the argument is true and create some kind of failure otherwise.
 */
namespace Checkers {
    
    /**
     * Checks the truthiness of the given boolean values through asserts.
     */
    struct Asserter {
        /**
         * Enforces truthiness of the given boolean value @b through asserts.
         *
         * @param b Boolean value whose truthiness should be enforced
         */
        static void enforce(bool b) {
            assert(b);
            (void) b; // prevent warnings in release-builds
        }
    };
    
    /**
     * Enforces truthiness of a given boolean value by throwing an exception in case of violation.
     */
    struct Enforcer {
        /**
         * Enforces that @a b is true and throws an Exception otherwise.
         *
         * @param b Boolean value whose truthiness should be enforced
         */
        static void enforce(bool b) {
            ::Aux::enforce(b);
        }
    };
    
    /**
     * Calls std::terminate if the bool is false.
     */
    struct Terminator {
        /**
         * Enforces truthiness of the given boolean value @b by terminating in case of violation.
         *
         * @param b Boolean value whose truthiness should be enforced
         */
        static void enforce(bool b) {
            if (!b) {
                std::terminate();
            }
        }
    };

    /**
     * Won't look at the bool (not even in debug-mode, which is how this differs from Asserter).
     */
    struct Ignorer {
        /**
         * Ignores the given boolean value and does not enforce anything.
         *
         * Useful for debugging purposes.
         *
         * @param Boolean whose value is ignored
         */
        static void enforce(bool) { }
    };
}

} // namespace Aux

#endif // NETWORKIT_AUXILIARY_ENFORCE_HPP_
