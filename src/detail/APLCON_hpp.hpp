#ifndef _APLCON_APLCON_HPP_HPP
#define _APLCON_APLCON_HPP_HPP 1

#include <type_traits>
#include <vector>
#include <functional>
#include <algorithm>

namespace APLCON_ {

// vectorize by template, that means wrap double value around initializer list {}

template<bool ReturnDouble>
struct vectorize_if {};

template<>
struct vectorize_if<true>  {
  static std::vector<double> get(const double& v) {
    return {v};
  }
};

template<>
struct vectorize_if<false>  {
  static std::vector<double> get(const std::vector<double>& v) {
    return v;
  }
};

// get some function traits like return type and number of arguments
// based on https://functionalcpp.wordpress.com/2013/08/05/function-traits/

template<template<typename,typename>class checker, typename... Ts>
struct is_all : std::true_type {};

template<template<typename,typename>class checker, typename T0, typename T1, typename... Ts>
struct is_all<checker, T0, T1, Ts...> :
    std::integral_constant< bool, checker<T0, T1>::value && is_all<checker, T0, Ts...>::value>
{};

// std::decay removes const and reference qualifiers
// which spoils the type comparison we actually want...
template<typename... Ts>
struct is_all_same_decayed : is_all<std::is_same, typename std::decay<Ts>::type ... > {};

// functor
template<class F>
struct function_traits
{
private:
  using call_type = function_traits<decltype(&F::operator())>; // how to prevent the compiler message here for non-callable types?

public:
  using return_type = typename call_type::return_type;

  static constexpr std::size_t arity = call_type::arity;
  static constexpr bool is_functor = true;

  template <typename Compare>
  struct all_args : call_type::template all_args<Compare> {};

};

template<class F>
struct function_traits<F&> : public function_traits<F> {};

template<class F>
struct function_traits<F&&> : public function_traits<F> {};

// anything else what is callable
template<typename R, typename... Args>
struct function_traits<R(Args...)>
{
  using return_type = R;

  static constexpr std::size_t arity = sizeof...(Args);
  static constexpr bool is_functor = false;

  template <typename Compare>
  struct all_args : is_all_same_decayed<Args..., Compare> {};

};

// function pointer
template<typename R, typename... Args>
struct function_traits<R(*)(Args...)> : public function_traits<R(Args...)> {};

// member function pointer
// ignore C, which makes all_args actually work...
template<class C, class R, class... Args>
struct function_traits<R(C::*)(Args...)> : public function_traits<R(Args...)> {};

// const member function pointer
// ignore C, which makes all_args actually work...
template<class C, class R, class... Args>
struct function_traits<R(C::*)(Args...) const> : public function_traits<R(Args...)> {};

// member object pointer
template<class C, class R>
struct function_traits<R(C::*)> : public function_traits<R(C&)> {};

// hide it from doxygen to prevent meaningless recursion warning
/// @cond
// this little template fun is called "pack of indices"
// it enables the nice definition of constraints via AddConstraint(...) method
// see http://stackoverflow.com/a/11044592
// and http://loungecpp.wikidot.com/tips-and-tricks%3aindices
template <std::size_t... Is>
struct indices {};

template <std::size_t N, std::size_t... Is>
struct build_indices : build_indices<N-1, N-1, Is...> {};

template <std::size_t... Is>
struct build_indices<0, Is...> : indices<Is...> {};
/// @endcond

} // end namespace APLCON_

#endif // _APLCON_APLCON_HPP_HPP

