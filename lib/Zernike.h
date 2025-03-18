#ifndef ZERNIKE_H
#define ZERNIKE_H

#include <cstdlib>  //std::size_t
#include <memory>   // unique_ptr

#include "Polynomial.h"
#include "util.h"

#ifndef MCT_MAX_ZERNIKE_ORDER
#define MCT_MAX_ZERNIKE_ORDER 40
#endif

template <int n, int m>
struct ZernikeRadial {
   private:
    template <int... l>
    static auto poly_even(std::integer_sequence<int, l...>) {
        return Polynomial<std::ratio<
            (n - l) % 2 == 0 && l >= m ? ((n - l) / 2 % 2 == 0 ? 1 : (-1)) *
                                             util::factorial(n - (n - l) / 2)
                                       : 0,
            (n - l) % 2 == 0
                ? util::factorial((n - l) / 2) *
                      util::factorial((n + m) / 2 - (n - l) / 2) *
                      util::factorial(util::abs((n - m) / 2 - (n - l) / 2))
                : 1>...>{};
    }

   public:
    using polynomial =
        std::conditional_t<(n - m) % 2 == 0,
                           decltype(poly_even(
                               std::make_integer_sequence<int, n + 1>{})),
                           Polynomial<std::ratio<0>>>;
};

namespace {

constexpr int index2to1(int n, int m) {
    return (n + 1) * (n + 1) / 4 + m / 2;
}

constexpr int index_n(int l) {
    if (l == 0) { return 0; }
    const int n = util::sqrt_int(4 * l) - 1;
    return index2to1(n, 0) + n / 2 < l ? n + 1 : n;
}

constexpr int index_m(int l) {
    const int n = index_n(l);
    return 2 * (l - index2to1(n, 0)) + n % 2;
}

struct ZernikeRadialWrapperBase {
    virtual double eval(double) const = 0;
    virtual ~ZernikeRadialWrapperBase() {}
};

template <int n, int m>
struct ZernikeRadialWrapper : public ZernikeRadialWrapperBase {
    double eval(double r) const override {
        return ZernikeRadial<n, m>::polynomial::eval(r);
    }
};

}  // namespace

static auto zernike_radial_at(int n, int m, double r) {
    constexpr auto radial_poly_count = index2to1(1 + MCT_MAX_ZERNIKE_ORDER, 0);
    static auto zernike_radials =
        ([]<auto... l>(std::integer_sequence<int, l...>)
             -> std::array<std::unique_ptr<ZernikeRadialWrapperBase>,
                           radial_poly_count> {
            return {std::make_unique<
                ZernikeRadialWrapper<index_n(l), index_m(l)>>()...};
        })(std::make_integer_sequence<int, radial_poly_count>{});
    return zernike_radials[static_cast<std::size_t>(index2to1(n, m))]->eval(r);
}
#endif  // ZERNIKE_H
