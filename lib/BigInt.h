#ifndef MCT_BIGINT_
#define MCT_BIGINT_

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
                  std::ratio_less<typename poly_leading<P1>::type,
                                  typename poly_leading<P2>::type>,
                  std::conjunction<
                      std::ratio_equal<typename poly_leading<P1>::type,
                                       typename poly_leading<P2>::type>,
                      bigint_less<typename poly_tail<P1>::type,
                                  typename poly_tail<P2>::type>>>>> {};

template <typename P1, typename P2>
struct bigint_less : bigint_less_helper<P1, P2> {};
template <typename C1, typename C2>
struct bigint_less<Polynomial<C1>, Polynomial<C2>> : std::ratio_less<C1, C2> {};
template <typename C>
struct bigint_less<Polynomial<C>, Polynomial<>>
    : std::ratio_less<C, std::ratio<0>> {};
template <typename C>
struct bigint_less<Polynomial<>, Polynomial<C>>
    : std::ratio_less<std::ratio<0>, C> {};

template <typename, typename>
struct bigint_equal;

template <typename P1, typename P2>
struct bigint_equal_helper
    : std::conjunction<std::integral_constant<bool, (P1::order == P2::order)>,
                       std::ratio_equal<typename poly_leading<P1>::type,
                                        typename poly_leading<P2>::type>,
                       bigint_equal<typename poly_tail<P1>::type,
                                    typename poly_tail<P2>::type>> {};

template <typename P1, typename P2>
struct bigint_equal : bigint_equal_helper<P1, P2> {};
template <typename C1, typename C2>
struct bigint_equal<Polynomial<C1>, Polynomial<C2>> : std::ratio_equal<C1, C2> {
};
template <typename C>
struct bigint_equal<Polynomial<C>, Polynomial<>>
    : std::ratio_equal<C, std::ratio<0>> {};
template <typename C>
struct bigint_equal<Polynomial<>, Polynomial<C>>
    : std::ratio_equal<std::ratio<0>, C> {};

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
template <typename C, typename P, bool to_positive = true>
struct normalize_bigint {};
template <bool to_positive, typename C0, typename C1, typename... Cs>
struct normalize_bigint<C0, Polynomial<C1, Cs...>, to_positive> {
    static constexpr std::intmax_t r = ((C0::num + C1::num) % bigint_base +
                                        (to_positive ? 1 : -1) * bigint_base) %
                                       bigint_base;
    using type = typename poly_prepend<
        std::ratio<r>,
        typename normalize_bigint<
            std::ratio<(C0::num + C1::num - r) / bigint_base>,
            Polynomial<Cs...>,
            to_positive>::type>::type;
};
template <bool to_positive, typename C>
struct normalize_bigint<C, Polynomial<>, to_positive> {
    static constexpr std::intmax_t r =
        (C::num % bigint_base + (to_positive ? 1 : -1) * bigint_base) %
        bigint_base;
    using type =
        Polynomial<std::ratio<r>, std::ratio<(C::num - r) / bigint_base>>;
};

template <typename P1, typename P2>
struct bigint_mul_impl {
    using type =
        typename normalize_bigint<std::ratio<0>, poly_mul<P1, P2>>::type;
};

template <typename C>
struct mixed_fraction_parts {
    static constexpr std::intmax_t integer_part = C::num / C::den;
    using fractional_part = std::ratio_subtract<C, std::ratio<integer_part>>;
};

template <typename P1, typename P2>
struct bigint_add_impl {
    using type =
        typename normalize_bigint<std::ratio<0>, poly_add<P1, P2>>::type;
};

template <typename P1, typename P2>
struct bigint_sub_impl {
    using type =
        typename normalize_bigint<std::ratio<0>,
                                  poly_sub<P1, P2>,
                                  bigint_less<P2, P1>::type::value>::type;
};

template <typename P>
struct collapse_leading_term {};
template <typename C, typename... Cs>
struct collapse_leading_term<Polynomial<C, Cs...>> {
    using tail_result = collapse_leading_term<Polynomial<Cs...>>;
    using carry = std::ratio<0>;

    using type =
        typename poly_prepend<std::ratio_add<C, typename tail_result::carry>,
                              typename tail_result::type>::type;
};
template <typename C>
struct collapse_leading_term<Polynomial<C>> {
    using type = Polynomial<>;
    using carry = std::ratio_multiply<C, std::ratio<bigint_base>>;
};

template <typename P1, typename P2>
struct bigint_div_helper {
    using type = std::conditional_t<
        std::ratio_greater_v<typename poly_leading<P1>::type,
                             typename poly_leading<P2>::type>,
        typename right_shift<
            Polynomial<std::ratio<poly_leading<P1>::type::num /
                                  (poly_leading<P2>::type::num + 1)>>,
            P1::order - P2::order>::type,
        typename bigint_div_helper<typename collapse_leading_term<P1>::type,
                                   P2>::type>;
};

template <typename P1, typename P2>
struct bigint_div_impl {
    using q = typename bigint_div_helper<P1, P2>::type;
    using sub_type = bigint_div_impl<
        typename bigint_sub_impl<P1,
                                 typename bigint_mul_impl<q, P2>::type>::type,
        P2>;

    using quotient =
        typename bigint_add_impl<q, typename sub_type::quotient>::type;
    using remainder = typename sub_type::remainder;
};

// template <typename P1, typename P2>
//     requires(bigint_less_equal_v<BigInt<P1>, BigInt<P2>>)
// struct bigint_div_impl<P1, P2> {};

template <std::intmax_t I>
struct to_bigint_impl {
    using type = typename poly_prepend<
        std::ratio<I % bigint_base>,
        typename to_bigint_impl<I / bigint_base>::type>::type;
};
template <>
struct to_bigint_impl<0> {
    using type = Polynomial<std::ratio<0>>;
};

template <typename C, typename P>
struct bigint_predecessor_helper {};
template <typename C0, typename C1, typename... Cs>
struct bigint_predecessor_helper<C0, Polynomial<C1, Cs...>> {
    using s = std::ratio_add<C1, C0>;
    using type = typename poly_prepend<
        std::ratio<(mixed_fraction_parts<s>::integer_part + bigint_base) %
                   bigint_base>,
        typename bigint_predecessor_helper<std::ratio<(s::num < 0 ? -1 : 0)>,
                                           Polynomial<Cs...>>::type>::type;
};
template <typename C>
struct bigint_predecessor_helper<C, Polynomial<>> {
    using type = Polynomial<C>;
};

template <typename P1, typename P2>
struct bigint_factorial_impl {
    using type = typename bigint_mul_impl<
        P1,
        typename bigint_factorial_impl<
            typename poly_trim<
                typename bigint_predecessor_helper<std::ratio<-1>,
                                                   P1>::type>::type,
            P2>::type>::type;
};
template <typename P>
struct bigint_factorial_impl<P, P> {
    using type = Polynomial<std::ratio<1>>;
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

// template <typename B1, typename B2>
// using bigint_div =
//     BigInt<typename impl::bigint_div_impl<typename B1::internal,
//                                           typename B2::internal>::type>;

template <std::intmax_t I>
using to_bigint = BigInt<
    typename impl::poly_trim<typename impl::to_bigint_impl<I>::type>::type>;

template <typename B1, typename B2>
using bigint_factorial_partial = BigInt<typename impl::poly_trim<
    typename impl::bigint_factorial_impl<typename B1::internal,
                                         typename B2::internal>::type>::type>;

template <typename B>
using bigint_factorial = bigint_factorial_partial<B, to_bigint<0>>;

#endif
