#ifndef ZERNIKE_H
#define ZERNIKE_H

#include <cstdlib>  //std::size_t
#include <memory>   // unique_ptr

#ifdef MCT_DEBUG_
#include <stdexcept>
#endif

#include "BigInt.h"
#include "Polynomial.h"
#include "util.h"

#ifndef MCT_MAX_ZERNIKE_POLAR_ORDER
#define MCT_MAX_ZERNIKE_POLAR_ORDER 20
#endif

namespace Zernike {

// OSA/ANSI standard indices
constexpr std::size_t index_l(auto n, auto m) {
    return static_cast<std::size_t>((n + 1) * n / 2 + (m + n) / 2);
}
constexpr int index_n(auto l) {
    return (util::sqrt_int(8 * static_cast<int>(l) + 1) - 1) / 2;
}
constexpr int index_m(auto l) {
    const int n = index_n(l);
    return 2 * static_cast<int>(l) - (n + 2) * n;
}
constexpr std::pair<int, int> index_nm(auto l) {
    const int n = index_n(l);
    return {n, 2 * static_cast<int>(l) - (n + 2) * n};
}

constexpr std::size_t basic_cap(auto mt) {
    return mt * (mt + 2);
}
constexpr std::size_t basic_index_l(auto n, auto m, auto mt) {
    return static_cast<std::size_t>(
        n <= mt ? index_l(n, m)
                : basic_cap(mt) - index_l(2 * mt - n, -static_cast<int>(m)));
}
constexpr int basic_index_n(auto l, auto mt) {
    const auto t = basic_cap(mt);
    return 2 * l < t ? inedx_n(l) : 2 * mt - index_n(t - l);
}
constexpr int basic_index_m(auto l, auto mt) {
    const auto t = basic_cap(mt);
    return 2 * l < t ? index_m(l) : -index_m(t - l);
}
constexpr std::pair<int, int> basic_index_nm(auto l, auto mt) {
    const auto t = basic_cap(mt);
    return 2 * static_cast<std::size_t>(l) < t
               ? index_nm(l)
               : (([mt](auto p) -> decltype(p) {
                     return {2 * mt - p.first, -p.second};
                 })(index_nm(t - l)));
}

namespace {

template <int n, int m>
struct RadialPolynomial {
   private:
    template <int... l>
    static auto poly_even(std::integer_sequence<int, l...>) {
        return typename impl::right_shift<
            Polynomial<std::ratio<
                (n - m - l) % 2 != 0
                    ? 0
                    : ((n - m - l) % 4 == 0 ? 1 : -1) *
                          from_bigint_v<bigint_div<
                              bigint_factorial_partial<
                                  to_bigint<(n + m + l) / 2>,
                                  to_bigint<(n - m - l) / 2>>,
                              bigint_mul<
                                  bigint_factorial<to_bigint<(l / 2 + m)>>,
                                  bigint_factorial<to_bigint<(l / 2)>>>>>>...>,
            m>::type{};
    }

   public:
    using polynomial =
        std::conditional_t<(n - m) % 2 == 0,
                           decltype(poly_even(
                               std::make_integer_sequence<int, n - m + 1>{})),
                           Polynomial<std::ratio<0>>>;
};

constexpr int radial_index_l(int n, int m) {
    return (n + 1) * (n + 1) / 4 + m / 2;
}

constexpr int radial_index_n(int l) {
    return util::sqrt_int(4 * l + 1) - 1;
}

constexpr int radial_index_m(int l) {
    const int n = radial_index_n(l);
    return 2 * l - n * (n + 2) / 2;
}

constexpr std::pair<int, int> radial_index_nm(int l) {
    const int n = radial_index_n(l);
    return {n, 2 * l - n * (n + 2) / 2};
}

// last index
constexpr int basic_radial_cap(int mt = MCT_MAX_ZERNIKE_POLAR_ORDER) {
    return mt * (mt + 3) / 2;
}

constexpr int basic_radial_index_l(int n,
                                   int m,
                                   int mt = MCT_MAX_ZERNIKE_POLAR_ORDER) {
    return n <= mt ? radial_index_l(n, m)
                   : basic_radial_cap(mt) - radial_index_l(2 * mt - n, m);
}

constexpr int basic_radial_index_n(int l,
                                   int mt = MCT_MAX_ZERNIKE_POLAR_ORDER) {
    const auto t = basic_radial_cap(mt);
    return 2 * l < t ? radial_index_n(l) : 2 * mt - radial_index_n(t - l);
}

constexpr int basic_radial_index_m(int l,
                                   int mt = MCT_MAX_ZERNIKE_POLAR_ORDER) {
    const auto t = basic_radial_cap(mt);
    return 2 * l < t ? radial_index_m(l) : radial_index_m(t - l);
}

constexpr std::pair<int, int> basic_radial_index_nm(
    int l,
    int mt = MCT_MAX_ZERNIKE_POLAR_ORDER) {
    const auto t = basic_radial_cap(mt);
    return 2 * l < t ? radial_index_nm(l) : (([mt](auto p) -> decltype(p) {
        return {2 * mt - p.first, p.second};
    })(radial_index_nm(t - l)));
}

struct WrapperBase {
    virtual double eval(double) const = 0;
    virtual double drvt(std::size_t, double) const = 0;
    virtual ~WrapperBase() {}
};

template <int n, int m>
struct Wrapper : public WrapperBase {
    double eval(double r) const override {
        return RadialPolynomial<n, m>::polynomial::eval(r);
    }
    double drvt(std::size_t d, double r) const override {
#ifdef MCT_DEBUG_
        if (d > 2) {
            throw std::runtime_error(
                "[Zernike::Wrapper] Derivative order can not exceed 2.");
        }
#endif
        return d == 0 ? eval(r)
               : d == 1
                   ? RadialPolynomial<n, m>::polynomial::template derivative<1>(
                         r)
                   : RadialPolynomial<n, m>::polynomial::template derivative<2>(
                         r);
    }
};

}  // namespace

double radial_at(int n, int m, double r, std::size_t d = 0);

template <typename T>
struct Series {
    using val_type = T;

    const std::size_t order;  // maximum polar order
    //  Maximum radial order = 2 * order

    Series(std::size_t polar_order = MCT_MAX_ZERNIKE_POLAR_ORDER)
        : order(polar_order), coef(basic_cap(polar_order) + 1) {}
    Series(std::size_t polar_order, std::vector<val_type> coefficients)
        : order(polar_order), coef{std::move(coefficients)} {
        const auto count = basic_cap(order) + 1;
        if (count != coef.size()) {
            std::cout << "[Zernike::Series] Specified order and the number of "
                         "provided coefficients do not match each order.\n";
            if (count < coef.size()) {
                std::cout << "[Zernike::Series] Coefficient list has "
                          << coef.size() << "elements and is truncted to "
                          << count << " elements.\n";
            } else {
                std::cout << "[Zernike::Series] Coefficient list has "
                          << coef.size()
                          << " elements and 0s are appended to make its size "
                          << count << ".\n";
            }
            coef.resize(count);
        }
    }

    template <util::Indexed2D U, util::Indexed1D V>
    Series(std::size_t polar_order,
           std::size_t radial_size,
           std::size_t polar_size,
           U data,
           const V& radial_coord)
        : Series(polar_order) {
        // force polar_order to be multiplier of 4 can reduce cache size to
        // polar_size/4
        const val_type theta_delta =
            2 * M_PI / static_cast<val_type>(polar_size);
        std::vector<std::array<double, 2>> sincos(polar_size);
        for (std::size_t i = 0; i < polar_size; ++i) {
            const auto theta = static_cast<val_type>(i) * theta_delta;
            sincos[i][0] = std::sin(theta);
            sincos[i][1] = std::cos(theta);
        }

        for (std::size_t l = 0; l < coef.size(); ++l) {
            const auto [n, m] = basic_index_nm(l, order);

            val_type r0 = 0.;
            val_type dr0 = 0;
            val_type dr1 = 0;
            double f0 = 0.;
            double f1 = 0.;
            double f2 = 0.;
            for (std::size_t j = 0; j < radial_size; ++j) {
                val_type circ_integral_l{};
                for (std::size_t i = 0; i < polar_size; ++i) {
                    const val_type simpson_coef =
                        ((i == 0 || i == polar_size - 1) && polar_size % 2 != 0
                             ? 5. / 6.
                             : static_cast<val_type>(2 * (1 + i % 2)) / 3.) *
                        theta_delta;
                    circ_integral_l +=
                        simpson_coef *
                        sincos[i * util::abs(m) % polar_size][m < 0 ? 0 : 1] *
                        data(j, i);
                }

                const val_type r = radial_coord[j];
                if (j % 2 == 0) {
                    dr0 = r - r0;
                    f1 = circ_integral_l * radial_at(n, util::abs(m), r) * r;
                } else {
                    dr1 = r - r0;
                    f2 = circ_integral_l * radial_at(n, util::abs(m), r) * r;
                    coef[l] += (dr0 + dr1) / 6. *
                               (2. * (f0 + f1 + f2) + dr0 / dr1 * (f1 - f2) +
                                dr1 / dr0 * (f1 - f0));
                    f0 = f2;
                }
                r0 = r;
            }

            coef[l] *= (n + 1) * (m == 0 ? 1. : 2.) / M_PI;
            for (std::size_t j = 0; j < radial_size; ++j) {
                const val_type r = radial_coord[j];
                for (std::size_t i = 0; i < polar_size; ++i) {
                    data(j, i) -=
                        coef[l] * radial_at(n, util::abs(m), r) *
                        sincos[i * util::abs(m) % polar_size][m < 0 ? 0 : 1];
                }
            }
        }
    }

    val_type operator()(val_type r, val_type theta) const {
        val_type f{};
        for (std::size_t l = 0; l < coef.size(); ++l) {
            const auto [n, m] = basic_index_nm(l, order);
            f += coef[l] * radial_at(n, m, r) *
                 (m == 0  ? 1.
                  : m < 0 ? std::sin(-m * theta)
                          : std::cos(m * theta));
        }
        return f;
    }

    val_type derivative(std::array<val_type, 2> pt,
                        std::array<std::size_t, 2> derivative_order) const {
        const auto [r, theta] = pt;
        const auto [rd, td] = derivative_order;

        val_type f{};
        for (std::size_t l = 0; l < coef.size(); ++l) {
            const auto [n, m] = basic_index_nm(l, order);
            const val_type t =
                util::abs(m) * theta + .5 * static_cast<val_type>(td) * M_PI;

            f += coef[l] * radial_at(n, m, r, rd) * std::pow(m, td) *
                 (m == 0  ? td > 0 ? 0. : 1.
                  : m < 0 ? std::sin(t)
                          : std::cos(t));
        }
        return f;
    }

    const auto& coefficient() const { return coef; }

   private:
    std::vector<val_type> coef;
};

}  // namespace Zernike
#endif  // ZERNIKE_H

#ifdef MCT_ZERNIKE_POLYNOMIAL_INSTANTIATION
namespace Zernike {
double radial_at(int n, int m, double r, std::size_t d) {
    constexpr auto radial_polynomial_count = 1 + basic_radial_cap();
    static auto zernike_radials =
        ([]<auto... l>(std::integer_sequence<int, l...>)
             -> std::array<std::unique_ptr<WrapperBase>,
                           radial_polynomial_count> {
            return {std::make_unique<Wrapper<basic_radial_index_n(l),
                                             basic_radial_index_m(l)>>()...};
        })(std::make_integer_sequence<int, radial_polynomial_count>{});
    return zernike_radials[static_cast<std::size_t>(basic_radial_index_l(n, m))]
        ->drvt(d, r);
}
}  // namespace Zernike
#endif
