#ifndef MEQ_POLYNOMIAL_H
#define MEQ_POLYNOMIAL_H

#include <algorithm>  //max
#include <array>
#include <cstdint>  // std::intmax_t
#include <cstdlib>  //std::size_t
#include <type_traits>

template <std::intmax_t... coefs>
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

template <std::intmax_t... cs, std::intmax_t... ps>
struct poly_derivative_helper<Polynomial<cs...>,
                              std::integer_sequence<std::intmax_t, ps...>> {
    using type = Polynomial<cs*(ps + 1)...>;
};

template <typename P, std::size_t n>
struct poly_derivative_impl {};
template <std::intmax_t c, std::intmax_t... cs>
struct poly_derivative_impl<Polynomial<c, cs...>, 0> {
    using type = Polynomial<c, cs...>;
};
template <std::size_t n>
struct poly_derivative_impl<Polynomial<>, n> {
    using type = Polynomial<0>;
};
template <std::size_t n, std::intmax_t c, std::intmax_t... cs>
struct poly_derivative_impl<Polynomial<c, cs...>, n> {
    using type = typename poly_derivative_impl<
        typename poly_derivative_helper<
            Polynomial<cs...>,
            std::make_integer_sequence<std::intmax_t, sizeof...(cs)>>::type,
        n - 1>::type;
};
};  // namespace impl

template <typename P, std::size_t n>
using poly_derivative = typename impl::poly_derivative_impl<P, n>::type;

template <std::intmax_t... coefs>
struct Polynomial {
    using Self = Polynomial<coefs...>;
    static constexpr std::size_t order = sizeof...(coefs);

    template <typename T>
    static constexpr T eval(T x) {
        return ([x]<std::intmax_t... cs>(Polynomial<cs...>) {
            return impl::eval_impl(std::array<T, order>{static_cast<T>(cs)...},
                                   x);
        })(Self{});
    }

    template <std::size_t n, typename T>
    static constexpr T derivative(T x) {
        return poly_derivative<Self, n>::eval(x);
    }

    template <typename OS>
    static void print_to(OS& out) {
        ([&out]<std::intmax_t... cs>(Polynomial<cs...>) {
            (..., (out << cs << ','));
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
template <std::intmax_t c2, std::intmax_t... c1s, std::intmax_t... c2s>
struct join<Polynomial<c1s...>, Polynomial<c2, c2s...>> {
    using type =
        typename join<Polynomial<c1s..., c2>, Polynomial<c2s...>>::type;
};

template <std::intmax_t c, std::size_t n, std::intmax_t... cs>
struct constant_list_impl : constant_list_impl<c, n - 1, c, cs...> {};

template <std::intmax_t c, std::intmax_t... cs>
struct constant_list_impl<c, 0, cs...> {
    using type = Polynomial<cs...>;
};

template <std::intmax_t c, std::size_t n>
struct constant_list {
    using type = typename constant_list_impl<c, n>::type;
};

template <typename P, std::size_t n>
struct pad {};

template <std::intmax_t c, std::intmax_t... cs, std::size_t n>
struct pad<Polynomial<c, cs...>, n> {
    using type =
        typename join<Polynomial<c>,
                      typename pad<Polynomial<cs...>, n - 1>::type>::type;
};
template <std::size_t n>
struct pad<Polynomial<>, n> {
    using type = typename constant_list<0, n>::type;
};

template <typename P1, typename P2>
struct poly_add_impl {};

template <std::intmax_t... c1s, std::intmax_t... c2s>
struct poly_add_impl<Polynomial<c1s...>, Polynomial<c2s...>> {
    using type = Polynomial<(c1s + c2s)...>;
};

template <typename P1, typename P2>
struct poly_sub_impl {};

template <std::intmax_t... c1s, std::intmax_t... c2s>
struct poly_sub_impl<Polynomial<c1s...>, Polynomial<c2s...>> {
    using type = Polynomial<(c1s - c2s)...>;
};

template <typename P, std::intmax_t c>
struct poly_scale_impl {};

template <std::intmax_t c, std::intmax_t... cs>
struct poly_scale_impl<Polynomial<cs...>, c> {
    using type = Polynomial<(cs * c)...>;
};

template <std::intmax_t c, typename P>
struct poly_prepend {};

template <std::intmax_t c, std::intmax_t... cs>
struct poly_prepend<c, Polynomial<cs...>> {
    using type = Polynomial<c, cs...>;
};

template <typename P>
struct poly_trim_helper {};

template <std::intmax_t c, std::intmax_t... cs>
struct poly_trim_helper<Polynomial<c, cs...>> {
    using tail_trim_result = poly_trim_helper<Polynomial<cs...>>;

    static constexpr bool trimming = tail_trim_result::trimming && c == 0;
    using type = std::conditional_t<
        trimming,
        Polynomial<>,
        typename poly_prepend<c, typename tail_trim_result::type>::type>;
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
                                    Polynomial<0>,
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

template <typename P, std::intmax_t c>
using poly_scale = typename impl::poly_scale_impl<P, c>::type;

namespace impl {

template <typename P, std::size_t n>
struct right_shift {
    using type = typename join<typename constant_list<0, n>::type, P>::type;
};

template <typename P1, typename P2>
struct poly_mul_impl {};

template <typename P2>
struct poly_mul_impl<Polynomial<>, P2> {
    using type = Polynomial<>;
};
template <std::intmax_t c1, std::intmax_t... c1s, typename P2>
struct poly_mul_impl<Polynomial<c1, c1s...>, P2> {
    using type =
        poly_add<poly_scale<P2, c1>,
                 typename right_shift<
                     typename poly_mul_impl<Polynomial<c1s...>, P2>::type,
                     1>::type>;
};

template <typename P>
struct poly_leading {};
template <>
struct poly_leading<Polynomial<>> {
    static constexpr std::intmax_t value = 0;
};
template <std::intmax_t c>
struct poly_leading<Polynomial<c>> {
    static constexpr std::intmax_t value = c;
};
template <std::intmax_t c, std::intmax_t... cs>
struct poly_leading<Polynomial<c, cs...>> {
    static constexpr std::intmax_t value =
        poly_leading<Polynomial<cs...>>::value;
};

template <typename P>
struct poly_tail {};
template <>
struct poly_tail<Polynomial<>> {
    using type = Polynomial<>;
};
template <std::intmax_t c>
struct poly_tail<Polynomial<c>> {
    using type = Polynomial<>;
};
template <std::intmax_t c, std::intmax_t... cs>
struct poly_tail<Polynomial<c, cs...>> {
    using type =
        typename poly_prepend<c, typename poly_tail<Polynomial<cs...>>::type>::
            type;
};

};  // namespace impl

template <typename P1, typename P2>
using poly_mul = typename impl::poly_mul_impl<P1, P2>::type;

#endif  // MEQ_POLYNOMIAL_H
