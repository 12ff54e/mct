#ifndef MEQ_BIGINT_H
#define MEQ_BIGINT_H

#include "Polynomial.h"

// Big int as polynomial
template <typename Poly>
struct BigInt : Poly {
    static constexpr std::intmax_t base = 1000;
    using internal = Poly;
};

static constexpr auto bigint_base = BigInt<Polynomial<>>::base;

namespace impl {

template <typename, typename>
struct bigint_less;

// Recursive alternatively helps short-circuiting the comparison
template <typename P1, typename P2>
struct bigint_less_helper
    : std::disjunction<
          std::integral_constant<bool, (P1::order < P2::order)>,
          std::conjunction<
              std::integral_constant<bool, (P1::order == P2::order)>,
              std::disjunction<
                  std::integral_constant<bool,
                                         (poly_leading<P1>::value <
                                          poly_leading<P2>::value)>,
                  std::conjunction<
                      std::integral_constant<bool,
                                             (poly_leading<P1>::value ==
                                              poly_leading<P2>::value)>,
                      bigint_less<typename poly_tail<P1>::type,
                                  typename poly_tail<P2>::type>>>>> {};

template <typename P1, typename P2>
struct bigint_less : bigint_less_helper<P1, P2> {};
template <std::intmax_t c1, std::intmax_t c2>
struct bigint_less<Polynomial<c1>, Polynomial<c2>>
    : std::integral_constant<bool, (c1 < c2)> {};
template <std::intmax_t c>
struct bigint_less<Polynomial<c>, Polynomial<>>
    : std::integral_constant<bool, (c < 0)> {};
template <std::intmax_t c>
struct bigint_less<Polynomial<>, Polynomial<c>>
    : std::integral_constant<bool, (0 < c)> {};

template <typename, typename>
struct bigint_equal;

template <typename P1, typename P2>
struct bigint_equal_helper
    : std::conjunction<std::integral_constant<bool, (P1::order == P2::order)>,
                       std::integral_constant<bool,
                                              (poly_leading<P1>::value ==
                                               poly_leading<P2>::value)>,
                       bigint_equal<typename poly_tail<P1>::type,
                                    typename poly_tail<P2>::type>> {};

template <typename P1, typename P2>
struct bigint_equal : bigint_equal_helper<P1, P2> {};
template <std::intmax_t c1, std::intmax_t c2>
struct bigint_equal<Polynomial<c1>, Polynomial<c2>>
    : std::integral_constant<bool, (c1 == c2)> {};
template <std::intmax_t c>
struct bigint_equal<Polynomial<c>, Polynomial<>>
    : std::integral_constant<bool, (c == 0)> {};
template <std::intmax_t c>
struct bigint_equal<Polynomial<>, Polynomial<c>>
    : std::integral_constant<bool, (c == 0)> {};

};  // namespace impl

template <typename B1, typename B2>
using bigint_less = typename impl::bigint_less<typename B1::internal,
                                               typename B2::internal>::type;
template <typename B1, typename B2>
constexpr bool bigint_less_v = bigint_less<B1, B2>::value;

template <typename B1, typename B2>
using bigint_equal = typename impl::bigint_equal<typename B1::internal,
                                                 typename B2::internal>::type;
template <typename B1, typename B2>
constexpr bool bigint_equal_v = bigint_equal<B1, B2>::value;

template <typename B1, typename B2>
using bigint_less_equal =
    std::disjunction<bigint_less<B1, B2>, bigint_equal<B1, B2>>;

template <typename B1, typename B2>
constexpr bool bigint_less_equal_v = bigint_less_equal<B1, B2>::value;

namespace impl {
template <std::intmax_t c, typename P, bool to_positive = true>
struct normalize_bigint {};
template <bool to_positive,
          std::intmax_t c0,
          std::intmax_t c1,
          std::intmax_t... cs>
struct normalize_bigint<c0, Polynomial<c1, cs...>, to_positive> {
    static constexpr std::intmax_t r =
        ((c0 + c1) % bigint_base + (to_positive ? 1 : -1) * bigint_base) %
        bigint_base;
    using type = typename poly_prepend<
        r,
        typename normalize_bigint<(c0 + c1 - r) / bigint_base,
                                  Polynomial<cs...>,
                                  to_positive>::type>::type;
};
template <bool to_positive, std::intmax_t c>
struct normalize_bigint<c, Polynomial<>, to_positive> {
    static constexpr std::intmax_t r =
        (c % bigint_base + (to_positive ? 1 : -1) * bigint_base) % bigint_base;
    using type = Polynomial<r, (c - r) / bigint_base>;
};

template <typename P1, typename P2>
struct bigint_mul_impl {
    using type = typename normalize_bigint<0, poly_mul<P1, P2>>::type;
};

template <typename P1, typename P2>
struct bigint_add_impl {
    using type = typename normalize_bigint<0, poly_add<P1, P2>>::type;
};

template <typename P1, typename P2>
struct bigint_sub_impl {
    using type =
        typename normalize_bigint<0,
                                  poly_sub<P1, P2>,
                                  bigint_less<P2, P1>::type::value>::type;
};

template <typename P>
struct collapse_leading_term {};
template <std::intmax_t c, std::intmax_t... cs>
struct collapse_leading_term<Polynomial<c, cs...>> {
    using tail_result = collapse_leading_term<Polynomial<cs...>>;
    static constexpr std::intmax_t carry = 0;

    using type = typename poly_prepend<(c + tail_result::carry),
                                       typename tail_result::type>::type;
};
template <std::intmax_t c>
struct collapse_leading_term<Polynomial<c>> {
    using type = Polynomial<>;
    static constexpr std::intmax_t carry = c * bigint_base;
};

template <typename P1, typename P2>
struct bigint_div_helper {
    using type = typename right_shift<
        Polynomial<(poly_leading<P1>::value / (poly_leading<P2>::value + 1))>,
        P1::order - P2::order>::type;
};

template <bool>
struct if_else {
    template <typename T, typename>
    using type = T;
};

template <>
struct if_else<false> {
    template <typename, typename T>
    using type = T;
};

template <typename P1, typename P2>
struct bigint_div_impl {
    using q =
        typename if_else<(poly_leading<P1>::value > poly_leading<P2>::value)>::
            template type<
                bigint_div_helper<P1, P2>,
                bigint_div_helper<typename collapse_leading_term<P1>::type,
                                  P2>>::type;

    using sub_type = bigint_div_impl<
        typename poly_trim<typename normalize_bigint<
            0,
            typename bigint_sub_impl<
                P1,
                typename bigint_mul_impl<q, P2>::type>::type>::type>::type,
        P2>;

    using quotient =
        typename bigint_add_impl<q, typename sub_type::quotient>::type;
    using remainder = typename sub_type::remainder;
};

template <typename P1, typename P2>
    requires(bigint_less_equal_v<BigInt<P1>, BigInt<P2>>)
struct bigint_div_impl<P1, P2> {
    using quotient = std::conditional_t<bigint_equal_v<BigInt<P1>, BigInt<P2>>,
                                        Polynomial<1>,
                                        Polynomial<0>>;
    using remainder = std::conditional_t<bigint_equal_v<BigInt<P1>, BigInt<P2>>,
                                         Polynomial<0>,
                                         P1>;
};

template <std::intmax_t I>
struct to_bigint_impl {
    using type = typename poly_prepend<
        I % bigint_base,
        typename to_bigint_impl<I / bigint_base>::type>::type;
};
template <>
struct to_bigint_impl<0> {
    using type = Polynomial<0>;
};

template <typename>
struct from_bigint_impl;
template <std::intmax_t c, std::intmax_t... cs>
struct from_bigint_impl<Polynomial<c, cs...>> {
    using type = std::integral_constant<
        std::intmax_t,
        c + bigint_base * from_bigint_impl<Polynomial<cs...>>::type::value>;
};
template <>
struct from_bigint_impl<Polynomial<>> {
    using type = std::integral_constant<std::intmax_t, 0>;
};

template <std::intmax_t c, typename P>
struct bigint_predecessor_helper {};
template <std::intmax_t c0, std::intmax_t c1, std::intmax_t... cs>
struct bigint_predecessor_helper<c0, Polynomial<c1, cs...>> {
    static constexpr std::intmax_t s = c0 + c1;
    using type = typename poly_prepend<
        (s + bigint_base) % bigint_base,
        typename bigint_predecessor_helper<(s < 0 ? -1 : 0),
                                           Polynomial<cs...>>::type>::type;
};
template <std::intmax_t c>
struct bigint_predecessor_helper<c, Polynomial<>> {
    using type = Polynomial<c>;
};

template <typename P1, typename P2>
struct bigint_factorial_impl {
    using type = typename bigint_mul_impl<
        P1,
        typename bigint_factorial_impl<
            typename poly_trim<
                typename bigint_predecessor_helper<-1, P1>::type>::type,
            P2>::type>::type;
};
template <typename P>
struct bigint_factorial_impl<P, P> {
    using type = Polynomial<1>;
};

};  // namespace impl

template <typename B1, typename B2>
using bigint_mul = BigInt<typename impl::poly_trim<
    typename impl::bigint_mul_impl<typename B1::internal,
                                   typename B2::internal>::type>::type>;

template <typename B1, typename B2>
using bigint_add = BigInt<typename impl::poly_trim<
    typename impl::bigint_add_impl<typename B1::internal,
                                   typename B2::internal>::type>::type>;

template <typename B1, typename B2>
using bigint_sub = BigInt<typename impl::poly_trim<
    typename impl::bigint_sub_impl<typename B1::internal,
                                   typename B2::internal>::type>::type>;

template <typename B1, typename B2>
using bigint_div = BigInt<typename impl::poly_trim<
    typename impl::bigint_div_impl<typename B1::internal,
                                   typename B2::internal>::quotient>::type>;

template <std::intmax_t I>
using to_bigint = BigInt<
    typename impl::poly_trim<typename impl::to_bigint_impl<I>::type>::type>;

template <typename B>
using from_bigint = typename impl::from_bigint_impl<typename B::internal>::type;

template <typename B>
constexpr std::intmax_t from_bigint_v = from_bigint<B>::value;

template <typename B1, typename B2>
using bigint_factorial_partial = BigInt<typename impl::poly_trim<
    typename impl::bigint_factorial_impl<typename B1::internal,
                                         typename B2::internal>::type>::type>;

template <typename B>
using bigint_factorial = bigint_factorial_partial<B, to_bigint<0>>;

#endif  // MEQ_BIGINT_H
