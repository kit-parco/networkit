// no-networkit-format
#ifndef NETWORKIT_AUXILIARY_FUNCTION_TRAITS_HPP_
#define NETWORKIT_AUXILIARY_FUNCTION_TRAITS_HPP_

#include <type_traits>
#include <tuple>

namespace Aux {

// Code taken from https://functionalcpp.wordpress.com/2013/08/05/function-traits/ and slightly modified

template<class F>
struct FunctionTraits;

// function pointer
template<class R, class... Args>
struct FunctionTraits<R(*)(Args...)> : public FunctionTraits<R(Args...)>
{};

template<class R, class... Args>
struct FunctionTraits<R(Args...)> {
    using result_type = R;

    static constexpr size_t arity = sizeof...(Args);

    template<size_t N, bool> struct arg_impl;
    template<size_t N> struct arg_impl<N, false> {};
    template<size_t N> struct arg_impl<N, true> {
        using type = typename std::tuple_element<N, std::tuple<Args...> >::type;
    };

    template <size_t N>
    using arg = arg_impl<N, N < arity>;
};


// member function pointer
template<class C, class R, class... Args>
struct FunctionTraits<R(C:: *)(Args...)> : public FunctionTraits<R(C &, Args...)>
{};

// const member function pointer
template<class C, class R, class... Args>
struct FunctionTraits<R(C:: *)(Args...) const> : public FunctionTraits<R(C &, Args...)>
{};

// member object pointer
template<class C, class R>
struct FunctionTraits<R(C:: *)> : public FunctionTraits<R(C &)>
{};

// functor
template<class F>
struct FunctionTraits {
private:
    using call_type = FunctionTraits<decltype(&F::operator())>;
public:
    using result_type = typename call_type::result_type;

    static constexpr size_t arity = call_type::arity - 1;

    template <size_t N>
    using arg = typename call_type::template arg < N + 1 >;
};

template<class F>
struct FunctionTraits<F &> : public FunctionTraits<F>
{};

template<class F>
struct FunctionTraits < F && > : public FunctionTraits<F>
{};



} // namespace Aux

#endif // NETWORKIT_AUXILIARY_FUNCTION_TRAITS_HPP_
