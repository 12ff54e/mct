#ifndef MEQ_ZERNIKE_H
#define MEQ_ZERNIKE_H

#include <cstdlib>  //std::size_t
#include <memory>   // unique_ptr

#ifdef MEQ_DEBUG_
#include <stdexcept>
#endif

#include "BigInt.h"
#include "Polynomial.h"
#include "util.h"

#ifndef MEQ_MAX_ZERNIKE_POLAR_ORDER
#define MEQ_MAX_ZERNIKE_POLAR_ORDER 20
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
    return static_cast<std::size_t>(mt * (mt + 2));
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
    const auto t = static_cast<decltype(l)>(basic_cap(mt));
    return 2 * l < t ? index_nm(l) : (([mt](auto p) -> decltype(p) {
        return {2 * static_cast<int>(mt) - p.first, -p.second};
    })(index_nm(t - l)));
}

namespace {

template <int n, int m>
struct RadialPolynomial {
   private:
    template <int... l>
    static auto poly_even(std::integer_sequence<int, l...>) {
        return typename impl::right_shift<
            Polynomial<
                (n - m - l) % 2 != 0
                    ? 0
                    : ((n - m - l) % 4 == 0 ? 1 : -1) *
                          from_bigint_v<bigint_div<
                              bigint_factorial_partial<
                                  to_bigint<(n + m + l) / 2>,
                                  to_bigint<(n - m - l) / 2>>,
                              bigint_mul<
                                  bigint_factorial<to_bigint<(l / 2 + m)>>,
                                  bigint_factorial<to_bigint<(l / 2)>>>>>...>,
            m>::type{};
    }

   public:
    using polynomial =
        std::conditional_t<(n - m) % 2 == 0,
                           decltype(poly_even(
                               std::make_integer_sequence<int, n - m + 1>{})),
                           Polynomial<0>>;
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
constexpr int basic_radial_cap(int mt = MEQ_MAX_ZERNIKE_POLAR_ORDER) {
    return mt * (mt + 3) / 2;
}

constexpr int basic_radial_index_l(int n,
                                   int m,
                                   int mt = MEQ_MAX_ZERNIKE_POLAR_ORDER) {
    return n <= mt ? radial_index_l(n, m)
                   : basic_radial_cap(mt) - radial_index_l(2 * mt - n, m);
}

constexpr int basic_radial_index_n(int l,
                                   int mt = MEQ_MAX_ZERNIKE_POLAR_ORDER) {
    const auto t = basic_radial_cap(mt);
    return 2 * l < t ? radial_index_n(l) : 2 * mt - radial_index_n(t - l);
}

constexpr int basic_radial_index_m(int l,
                                   int mt = MEQ_MAX_ZERNIKE_POLAR_ORDER) {
    const auto t = basic_radial_cap(mt);
    return 2 * l < t ? radial_index_m(l) : radial_index_m(t - l);
}

constexpr std::pair<int, int> basic_radial_index_nm(
    int l,
    int mt = MEQ_MAX_ZERNIKE_POLAR_ORDER) {
    const auto t = basic_radial_cap(mt);
    return 2 * l < t ? radial_index_nm(l) : (([mt](auto p) -> decltype(p) {
        return {2 * mt - p.first, p.second};
    })(radial_index_nm(t - l)));
}

struct WrapperBase {
    virtual double eval(double) const = 0;
    virtual ~WrapperBase() {}
};

template <int n, int m>
struct Wrapper : public WrapperBase {
    double eval(double r) const override {
        return RadialPolynomial<n, m>::polynomial::eval(r);
    }
};

}  // namespace

double radial_at(int n, int m, double r);

template <typename T>
struct Series {
    using val_type = T;

    const int order;  // maximum polar order
    //  Maximum radial order = 2 * order

    Series(int polar_order = MEQ_MAX_ZERNIKE_POLAR_ORDER)
        : order(polar_order), coef(basic_cap(polar_order) + 1) {}
    Series(int polar_order, std::vector<val_type> coefficients)
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
    Series(int polar_order,
           std::size_t radial_size,
           std::size_t polar_size,
           const U& data,
           const V& radial_coord)
        : Series(polar_order) {
        // force polar_size to be multiplier of 4 can reduce cache size to
        // polar_size/4
        const val_type theta_delta =
            2 * M_PI / static_cast<val_type>(polar_size);
        std::vector<std::array<val_type, 2>> sincos(polar_size);
        for (std::size_t i = 0; i < polar_size; ++i) {
            const auto theta = static_cast<val_type>(i) * theta_delta;
            sincos[i][0] = std::sin(theta);
            sincos[i][1] = std::cos(theta);
        }

        // stores \int_0^{2\pi}d\theta f(r,\theta)*trig(m*\theta), for each m
        // and r
        std::vector<val_type> angular_integrals(
            radial_size * static_cast<std::size_t>(2 * order + 1));
        const auto calc_angular_integral = [&](auto n, auto m, auto ri,
                                               auto r) -> val_type {
            auto& val =
                angular_integrals[radial_size *
                                      static_cast<std::size_t>(m + order) +
                                  ri];
            if (n == util::abs(m)) {
                for (std::size_t i = 0; i < polar_size; ++i) {
                    const auto angular_val =
                        data(ri, i) *
                        sincos[i * static_cast<std::size_t>(util::abs(m)) %
                               polar_size][m < 0 ? 0 : 1];
                    const val_type simpson_coef =
                        ((i == 0 || i == polar_size - 1) && polar_size % 2 != 0
                             ? 5. / 6.
                             : static_cast<val_type>(2 * (1 + i % 2)) / 3.) *
                        theta_delta;
                    val += simpson_coef * angular_val;
                }
            } else {
                val -= coef[basic_index_l(n - 2, m, order)] *
                       radial_at(n - 2, util::abs(m), r) * M_PI *
                       (m == 0 ? 2. : 1.);
            }
            return val;
        };

        for (std::size_t l = 0; l < coef.size(); ++l) {
            const auto [n, m] = basic_index_nm(l, order);

            val_type r0 = 0.;
            val_type f0 = 0.;

            std::array<val_type, 2> dr;
            std::array<val_type, 2> fr;
            for (std::size_t j = 0; j < radial_size; ++j) {
                const val_type r = radial_coord[j];
                const auto radial_val = radial_at(n, util::abs(m), r);
                const auto angular_integral = calc_angular_integral(n, m, j, r);

                dr[j % 2] = r - r0;
                fr[j % 2] = angular_integral * radial_val * r;
                if (j % 2 != 0) {
                    coef[l] += (dr[0] + dr[1]) / 6. *
                               (2. * (f0 + fr[0] + fr[1]) +
                                dr[0] / dr[1] * (fr[0] - fr[1]) +
                                dr[1] / dr[0] * (fr[0] - f0));
                    f0 = fr[1];
                }
                r0 = r;
            }

            coef[l] *= (n + 1) * (m == 0 ? 1. : 2.) / M_PI;
        }
    }

    val_type operator()(val_type r, val_type theta) const {
        val_type f{};
        std::array<val_type, 2 * MEQ_MAX_ZERNIKE_POLAR_ORDER + 1> trig_buffer;
        for (int i = 0; i <= 2 * order; ++i) {
            trig_buffer[static_cast<std::size_t>(i)] =
                i < order ? std::sin((order - i) * theta)
                          : std::cos((i - order) * theta);
        }
        for (std::size_t l = 0; l < coef.size(); ++l) {
            const auto [n, m] = basic_index_nm(l, order);
            f += coef[l] * radial_at(n, util::abs(m), r) *
                 trig_buffer[static_cast<std::size_t>(m + order)];
        }
        return f;
    }

    const auto& coefficient() const { return coef; }

   private:
    std::vector<val_type> coef;
};

}  // namespace Zernike
#endif  // ZERNIKE_H

#ifdef MEQ_ZERNIKE_POLYNOMIAL_INSTANTIATION
namespace Zernike {
double radial_at(int n, int m, double r) {
    constexpr auto radial_polynomial_count = 1 + basic_radial_cap();
    static auto zernike_radials =
        ([]<auto... l>(std::integer_sequence<int, l...>)
             -> std::array<std::unique_ptr<WrapperBase>,
                           radial_polynomial_count> {
            return {std::make_unique<Wrapper<basic_radial_index_n(l),
                                             basic_radial_index_m(l)>>()...};
        })(std::make_integer_sequence<int, radial_polynomial_count>{});
    return zernike_radials[static_cast<std::size_t>(basic_radial_index_l(n, m))]
        ->eval(r);
}
}  // namespace Zernike
#endif  // MEQ_ZERNIKE_H
