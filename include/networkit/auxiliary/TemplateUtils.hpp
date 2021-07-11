// no-networkit-format
#ifndef NETWORKIT_AUXILIARY_TEMPLATE_UTILS_HPP_
#define NETWORKIT_AUXILIARY_TEMPLATE_UTILS_HPP_

#include <type_traits>

/**
 * The macro makes the use of std::enable_if much easier by removing all the boilerplate.
 *
 * The 'what' argument must be a valid identifier that describes, what is required.
 * Example: you want to require that T is an integer, than what should be something like
 *          T_is_integer
 *
 * The second argument must be a constant boolean expression. In order to allow templates that
 * would otherwise confuse the preprocessor, this is passed as variadic argument (this is however
 * an implementation detail and you should never actually pass more than one argument.
 */
#define AUX_REQUIRE(what, ...) class what = typename ::std::enable_if<__VA_ARGS__>::type


/**
 * If two templates only differ in the default-values of their template-arguments
 * C++ considers this to be a compilation-error. A simple way to prevent this is to
 * add further defaulted template-arguments, which is what this macro provides in a
 * semantic way.
 *
 * Example:
 *
 * template<typename T, AUX_REQUIRE(isInteger, [...])>
 * void fun(T) {...}
 *
 * template<typename T, AUX_REQUIRE(isFloat, [...]), AUX_DUMMY_ARGUMENT>
 * void fun(T) {...}
 *
 * Note however that it may often be a cleaner solution to use compile-time dispatching instead
 * of hacks that envolve enable-if.
 */
#define AUX_DUMMY_ARGUMENT class = void

namespace Aux {

/**
 * This is a backport of C++14 std::decay_t
 */
template<typename T> using decay_t = typename std::decay<T>::type;

/**
 * Returns the corresponding std::integral_constant<bool, B> to a boolean
 * value B.
 */
template<bool B> using boolToType = std::integral_constant<bool, B>;

/**
 * This is just a wrapper around std::is_same in order to provide a nicer interface.
 *
 * With C++14 this could use template-variables, but with C++11 we are stuck with constexpr-functions.
 */
template<typename T1, typename T2>
constexpr bool isSame() {
    return std::is_same<T1, T2>::value;
}

/**
 * Checks whether Derived is either a type that derives from Base or is identical to Base.
 *
 * This is very usefull for situations in which you work with iterator-traits or the like.
 */
template<typename Base, typename Derived>
constexpr bool isBaseOrSame() {
    return isSame<Base, Derived>() || std::is_base_of<Base, Derived>::value;
}

} // namespace Aux

#endif // NETWORKIT_AUXILIARY_TEMPLATE_UTILS_HPP_
