#ifndef MCT_BIGINT_
#define MCT_BIGINT_

#include "Polynomial.h"

// Big int as polynomial
template <typename Poly>
struct BigInt : Poly {
    static constexpr std::intmax_t base = 1000;
    using internal = Poly;
};

namespace impl {

template <typename C, typename P>
struct normalize_bigint {};
template <typename C0, typename C1, typename... Cs>
struct normalize_bigint<C0, Polynomial<C1, Cs...>> {
    using type = typename poly_prepend<
        std::ratio<(C0::num + C1::num) % BigInt<Polynomial<>>::base>,
        typename normalize_bigint<
            std::ratio<(C0::num + C1::num) / BigInt<Polynomial<>>::base>,
            Polynomial<Cs...>>::type>::type;
};
template <typename C>
struct normalize_bigint<C, Polynomial<>> {
    using type = Polynomial<std::ratio<C::num % BigInt<Polynomial<>>::base>,
                            std::ratio<C::num / BigInt<Polynomial<>>::base>>;
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

template <typename P>
struct normalize_bigint_frac {};
template <typename C0, typename... Cs>
struct normalize_bigint_frac<Polynomial<C0, Cs...>> {
    using tail_result = normalize_bigint_frac<Polynomial<Cs...>>;
    using sum = mixed_fraction_parts<std::ratio_add<
        C0,
        std::ratio_multiply<typename tail_result::fraction,
                            std::ratio<BigInt<Polynomial<>>::base>>>>;
    using fraction = typename sum::fractional_part;
    using type = typename poly_prepend<std::ratio<sum::integer_part>,
                                       typename tail_result::type>::type;
};
template <>
struct normalize_bigint_frac<Polynomial<>> {
    using fraction = std::ratio<0>;
    using type = Polynomial<>;
};

template <typename P1, typename P2>
struct bigint_div_impl {
    using poly_divmod_result = poly_divmod<P1, P2>;
    using q = normalize_bigint_frac<typename poly_divmod_result::quotient>;
    using r = normalize_bigint_frac<typename poly_divmod_result::remainder>;

    using type = poly_add<
        typename q::type,
        poly_div<
            typename normalize_bigint<
                std::ratio<mixed_fraction_parts<std::ratio_add<
                    std::ratio_multiply<typename q::fraction,
                                        std::ratio<BigInt<Polynomial<>>::base>>,
                    typename r::fraction>>::integer_part>,
                typename r::type>::type,
            P2>>;
};

template <std::intmax_t I>
struct to_bigint_impl {
    using type = typename poly_prepend<
        std::ratio<I % BigInt<Polynomial<>>::base>,
        typename to_bigint_impl<I / BigInt<Polynomial<>>::base>::type>::type;
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
        std::ratio<(mixed_fraction_parts<s>::integer_part +
                    BigInt<Polynomial<>>::base) %
                   BigInt<Polynomial<>>::base>,
        typename bigint_predecessor_helper<std::ratio<(s::num < 0 ? -1 : 0)>,
                                           Polynomial<Cs...>>::type>::type;
};
template <typename C>
struct bigint_predecessor_helper<C, Polynomial<>> {
    using type = Polynomial<C>;
};

template <typename P>
struct bigint_factorial_impl {
    using type = typename bigint_mul_impl<
        P,
        typename bigint_factorial_impl<typename poly_trim<
            typename bigint_predecessor_helper<std::ratio<-1>,
                                               P>::type>::type>::type>::type;
};
template <>
struct bigint_factorial_impl<Polynomial<std::ratio<0>>> {
    using type = Polynomial<std::ratio<1>>;
};

};  // namespace impl

template <typename B1, typename B2>
using bigint_mul = BigInt<typename impl::poly_trim<
    typename impl::bigint_mul_impl<typename B1::internal,
                                   typename B2::internal>::type>::type>;

template <typename B1, typename B2>
using bigint_div =
    BigInt<typename impl::bigint_div_impl<typename B1::internal,
                                          typename B2::internal>::type>;

template <std::intmax_t I>
using to_bigint = BigInt<
    typename impl::poly_trim<typename impl::to_bigint_impl<I>::type>::type>;

template <typename B>
using bigint_factorial = BigInt<typename impl::poly_trim<
    typename impl::bigint_factorial_impl<typename B::internal>::type>::type>;

#endif
