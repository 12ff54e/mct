#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <algorithm>  //max
#include <array>
#include <cstdlib>  //std::size_t
#include <ratio>
#include <type_traits>

namespace {

template <typename A, typename T, std::size_t... idx>
constexpr T eval_impl_odd(A, T, std::index_sequence<idx...>);
template <typename A, typename T, std::size_t... idx>
constexpr T eval_impl_even(A, T, std::index_sequence<idx...>);

template <typename A, typename T>
constexpr T eval_impl(A arr, T pows) {
    return arr.size() == 1 ? arr[0]
           : arr.size() % 2 == 0
               ? eval_impl_even(arr, pows,
                                std::make_index_sequence<arr.size() / 2>{})
               : eval_impl_odd(arr, pows,
                               std::make_index_sequence<arr.size() / 2>{});
}

template <typename A, typename T, std::size_t... idx>
constexpr T eval_impl_odd(A arr, T pows, std::index_sequence<idx...>) {
    return eval_impl(
        std::array<T, sizeof...(idx) + 1>{
            (arr[2 * idx] + arr[2 * idx + 1] * pows)..., arr[arr.size() - 1]},
        pows * pows);
}
template <typename A, typename T, std::size_t... idx>
constexpr T eval_impl_even(A arr, T pows, std::index_sequence<idx...>) {
    return eval_impl(std::array<T, sizeof...(idx)>{(
                         arr[2 * idx] + arr[2 * idx + 1] * pows)...},
                     pows * pows);
}

};  // namespace

template <typename... Coefs>
struct Polynomial {
    using Self = Polynomial<Coefs...>;
    static constexpr std::size_t order = sizeof...(Coefs);

    template <typename T>
    static constexpr T eval(T x) {
        return ([x]<typename... Cs>(Polynomial<Cs...>) {
            return ::eval_impl(
                std::array<T, order>{
                    (static_cast<T>(Cs::num) / static_cast<T>(Cs::den))...},
                x);
        })(Self{});
    }

    template <typename OS>
    void print_to(OS& out) {
        ([&out]<typename... Cs>(Polynomial<Cs...>) {
            (..., (out << Cs::num << "/" << Cs::den << ","));
        })(Self{});
    }
};

namespace {

template <typename P1, typename P2>
struct join {};

template <typename P1>
struct join<P1, Polynomial<>> {
    using type = P1;
};
template <typename C2, typename... C1s, typename... C2s>
struct join<Polynomial<C1s...>, Polynomial<C2, C2s...>> {
    using type =
        typename join<Polynomial<C1s..., C2>, Polynomial<C2s...>>::type;
};

template <typename T, std::size_t n, typename... Ts>
struct constant_list_impl : constant_list_impl<T, n - 1, T, Ts...> {};

template <typename T, typename... Ts>
struct constant_list_impl<T, 0, Ts...> {
    using type = Polynomial<Ts...>;
};

template <typename T, std::size_t n>
struct constant_list {
    using type = typename constant_list_impl<T, n>::type;
};

template <typename P, std::size_t n>
struct pad {};

template <typename C, typename... Cs, std::size_t n>
struct pad<Polynomial<C, Cs...>, n> {
    using type =
        typename join<Polynomial<C>,
                      typename pad<Polynomial<Cs...>, n - 1>::type>::type;
};
template <std::size_t n>
struct pad<Polynomial<>, n> {
    using type = typename constant_list<std::ratio<0>, n>::type;
};

template <typename P, std::size_t n>
struct right_shift {
    using type =
        typename join<typename constant_list<std::ratio<0>, n>::type, P>::type;
};

template <typename P1, typename P2>
struct poly_add_impl {};

template <typename... C1s, typename... C2s>
struct poly_add_impl<Polynomial<C1s...>, Polynomial<C2s...>> {
    using type = Polynomial<std::ratio_add<C1s, C2s>...>;
};

template <typename P1, typename P2>
struct poly_sub_impl {};

template <typename... C1s, typename... C2s>
struct poly_sub_impl<Polynomial<C1s...>, Polynomial<C2s...>> {
    using type = Polynomial<std::ratio_subtract<C1s, C2s>...>;
};

template <typename P, typename T>
struct poly_scale_impl {};

template <typename T, typename... Cs>
struct poly_scale_impl<Polynomial<Cs...>, T> {
    using type = Polynomial<std::ratio_multiply<Cs, T>...>;
};

};  // namespace

template <typename P1, typename P2>
using poly_add = typename ::poly_add_impl<
    typename ::pad<P1, std::max(P1::order, P2::order)>::type,
    typename ::pad<P2, std::max(P1::order, P2::order)>::type>::type;

template <typename P1, typename P2>
using poly_sub = typename ::poly_sub_impl<
    typename ::pad<P1, std::max(P1::order, P2::order)>::type,
    typename ::pad<P2, std::max(P1::order, P2::order)>::type>::type;

template <typename P, typename T>
using poly_scale = typename ::poly_scale_impl<P, T>::type;

namespace {

template <typename P1, typename P2>
struct poly_mul_impl {};

template <typename P2>
struct poly_mul_impl<Polynomial<>, P2> {
    using type = Polynomial<>;
};
template <typename C1, typename... C1s, typename P2>
struct poly_mul_impl<Polynomial<C1, C1s...>, P2> {
    using type =
        poly_add<poly_scale<P2, C1>,
                 typename right_shift<
                     typename poly_mul_impl<Polynomial<C1s...>, P2>::type,
                     1>::type>;
};

};  // namespace

template <typename P1, typename P2>
using poly_mul = typename ::poly_mul_impl<P1, P2>::type;

#endif  // POLYNOMIAL_H
