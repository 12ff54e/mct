#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <algorithm>  //max
#include <array>
#include <cstdlib>  //std::size_t
#include <ratio>
#include <type_traits>

template <typename... Coefs>
struct Polynomial;

namespace impl {

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
    static_cast<void>(arr);  // avoid unused-but-set-parameter warning in gcc
    return eval_impl(std::array<T, sizeof...(idx)>{(
                         arr[2 * idx] + arr[2 * idx + 1] * pows)...},
                     pows * pows);
}

template <typename P, typename IS>
struct poly_derivative_helper {};

template <typename... Cs, std::intmax_t... ps>
struct poly_derivative_helper<Polynomial<Cs...>,
                              std::integer_sequence<std::intmax_t, ps...>> {
    using type = Polynomial<std::ratio_multiply<Cs, std::ratio<ps + 1>>...>;
};

template <typename P, std::size_t n>
struct poly_derivative_impl {};
template <typename C, typename... Cs>
struct poly_derivative_impl<Polynomial<C, Cs...>, 0> {
    using type = Polynomial<C, Cs...>;
};
template <std::size_t n>
struct poly_derivative_impl<Polynomial<>, n> {
    using type = Polynomial<std::ratio<0, 1>>;
};
template <std::size_t n, typename C, typename... Cs>
struct poly_derivative_impl<Polynomial<C, Cs...>, n> {
    using type = typename poly_derivative_impl<
        typename poly_derivative_helper<
            Polynomial<Cs...>,
            std::make_integer_sequence<std::intmax_t, sizeof...(Cs)>>::type,
        n - 1>::type;
};
};  // namespace impl

template <typename P, std::size_t n>
using poly_derivative = typename impl::poly_derivative_impl<P, n>::type;

template <typename... Coefs>
struct Polynomial {
    using Self = Polynomial<Coefs...>;
    static constexpr std::size_t order = sizeof...(Coefs);

    template <typename T>
    static constexpr T eval(T x) {
        return ([x]<typename... Cs>(Polynomial<Cs...>) {
            return impl::eval_impl(
                std::array<T, order>{
                    (static_cast<T>(Cs::num) / static_cast<T>(Cs::den))...},
                x);
        })(Self{});
    }

    template <std::size_t n, typename T>
    static constexpr T derivative(T x) {
        return poly_derivative<Self, n>::eval(x);
    }

    template <typename OS>
    static void print_to(OS& out) {
        ([&out]<typename... Cs>(Polynomial<Cs...>) {
            const auto print_ratio = [&out]<typename C>(C) {
                out << C::num;
                if (C::den != 1) { out << '/' << C::den; }
                out << ',';
            };
            (..., print_ratio(Cs{}));
        })(Self{});
    }
};

namespace impl {

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

template <typename C, typename P>
struct poly_prepend {};

template <typename C, typename... Cs>
struct poly_prepend<C, Polynomial<Cs...>> {
    using type = Polynomial<C, Cs...>;
};

template <typename P>
struct poly_trim_helper {};

template <typename C, typename... Cs>
struct poly_trim_helper<Polynomial<C, Cs...>> {
    using tail_trim_result = poly_trim_helper<Polynomial<Cs...>>;

    static constexpr bool trimming =
        tail_trim_result::trimming && std::is_same_v<C, std::ratio<0>>;
    using type = std::conditional_t<
        trimming,
        Polynomial<>,
        typename poly_prepend<C, typename tail_trim_result::type>::type>;
};

template <>
struct poly_trim_helper<Polynomial<>> {
    static constexpr bool trimming = true;
    using type = Polynomial<>;
};

template <typename P>
struct poly_trim {
    using helper_result = typename poly_trim_helper<P>::type;
    using type = std::conditional_t<std::is_same_v<helper_result, Polynomial<>>,
                                    Polynomial<std::ratio<0>>,
                                    helper_result>;
};

};  // namespace impl

template <typename P1, typename P2>
using poly_add = typename impl::poly_trim<typename impl::poly_add_impl<
    typename impl::pad<P1, std::max(P1::order, P2::order)>::type,
    typename impl::pad<P2, std::max(P1::order, P2::order)>::type>::type>::type;

template <typename P1, typename P2>
using poly_sub = typename impl::poly_trim<typename impl::poly_sub_impl<
    typename impl::pad<P1, std::max(P1::order, P2::order)>::type,
    typename impl::pad<P2, std::max(P1::order, P2::order)>::type>::type>::type;

template <typename P, typename T>
using poly_scale = typename impl::poly_scale_impl<P, T>::type;

namespace impl {

template <typename P, std::size_t n>
struct right_shift {
    using type =
        typename join<typename constant_list<std::ratio<0>, n>::type, P>::type;
};

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

template <typename P>
struct poly_leading {};

template <>
struct poly_leading<Polynomial<>> {
    using type = std::ratio<0>;
};
template <typename C>
struct poly_leading<Polynomial<C>> {
    using type = C;
};
template <typename C, typename... Cs>
struct poly_leading<Polynomial<C, Cs...>> {
    using type = typename poly_leading<Polynomial<Cs...>>::type;
};

template <typename P1, typename P2>
struct poly_div_helper {
    using type = typename right_shift<
        Polynomial<std::ratio_divide<typename poly_leading<P1>::type,
                                     typename poly_leading<P2>::type>>,
        P1::order - P2::order>::type;
};

template <typename P1, typename P2>
struct poly_div_impl {
    using q = typename poly_div_helper<P1, P2>::type;
    using tail_result =
        poly_div_impl<poly_sub<P1, typename poly_mul_impl<q, P2>::type>, P2>;

    using quotient = poly_add<q, typename tail_result::quotient>;
    using remainder = typename tail_result::remainder;
};
template <typename P1, typename P2>
    requires(P1::order < P2::order ||
             std::is_same_v<P1, Polynomial<std::ratio<0>>>)
struct poly_div_impl<P1, P2> {
    using quotient = Polynomial<std::ratio<0>>;
    using remainder = P1;
};

};  // namespace impl

template <typename P1, typename P2>
using poly_mul = typename impl::poly_mul_impl<P1, P2>::type;

template <typename P1, typename P2>
using poly_div = typename impl::poly_div_impl<P1, P2>::quotient;

template <typename P1, typename P2>
using poly_mod = typename impl::poly_div_impl<P1, P2>::remainder;

template <typename P1, typename P2>
using poly_divmod = impl::poly_div_impl<P1, P2>;

#endif  // POLYNOMIAL_H
