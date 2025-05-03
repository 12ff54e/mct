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

#ifndef MCT_MAX_ZERNIKE_ORDER
#define MCT_MAX_ZERNIKE_ORDER 14
#endif

namespace Zernike {

// OSA/ANSI standard indices
constexpr int index_l(int n, int m) {
    return (n + 1) * n / 2 + (m + n) / 2;
}
constexpr int index_n(int l) {
    return (util::sqrt_int(8 * l + 1) - 1) / 2;
}
constexpr int index_m(int l) {
    const int n = index_n(l);
    return 2 * l - (n + 2) * n;
}
constexpr std::pair<int, int> index_nm(int l) {
    const int n = index_n(l);
    return {n, 2 * l - (n + 2) * n};
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

constexpr int radial_index2to1(int n, int m) {
    return (n + 1) * (n + 1) / 4 + m / 2;
}

constexpr int radial_index_n(int l) {
    return util::sqrt_int(4 * l + 1) - 1;
}

constexpr int radial_index_m(int l) {
    const int n = radial_index_n(l);
    return 2 * (l - radial_index2to1(n, 0)) + n % 2;
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

    const std::size_t order;

    Series(std::size_t radial_order = MCT_MAX_ZERNIKE_ORDER)
        : order(radial_order) {}
    Series(std::size_t radial_order, std::vector<val_type> coefficients)
        : order(radial_order), coef{std::move(coefficients)} {
        const int count = index_l(radial_order, radial_order) + 1;
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

    val_type operator()(val_type r, val_type theta) const {
        val_type f{};
        for (std::size_t l = 0; l < coef.size(); ++l) {
            const auto [n, m] = index_nm(static_cast<int>(l));
            f += coef[l] * radial_at(n, m, r) *
                 (m == 0  ? 1.
                  : m < 0 ? std::sin(m * theta)
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
            const auto [n, m] = index_nm(static_cast<int>(l));
            const double t = m * theta + .5 * td * M_PI;

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
    constexpr auto radial_polynomial_count =
        radial_index2to1(1 + MCT_MAX_ZERNIKE_ORDER, 0);
    static auto zernike_radials =
        ([]<auto... l>(std::integer_sequence<int, l...>)
             -> std::array<std::unique_ptr<WrapperBase>,
                           radial_polynomial_count> {
            return {std::make_unique<
                Wrapper<radial_index_n(l), radial_index_m(l)>>()...};
        })(std::make_integer_sequence<int, radial_polynomial_count>{});
    return zernike_radials[static_cast<std::size_t>(radial_index2to1(n, m))]
        ->drvt(d, r);
}
}  // namespace Zernike
#endif
