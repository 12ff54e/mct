#include <cmath>
#include <iostream>
#include <limits>

#include "Assertion.hpp"
#include "util.h"

/**
 * @brief Find out if two float point number is approximately equal (diff up to
 * 3 last digits)
 */
template <typename T>
bool approx_eq(T a, T b) {
    return std::abs(a - b) <= T{4} * std::numeric_limits<T>::epsilon() *
                                  std::min(std::abs(a), std::abs(b));
}

int main() {
    Assertion assertion;

    long double sqrt2 =
        util::find_root([](long double x) { return x * x - 2; }, 1.l, 2.l);
    assertion(approx_eq(sqrt2, M_SQRT2l),
              std::string{"Find root function not correct! Error ~ "} +
                  std::to_string(sqrt2 - M_SQRT2l));

    double integral =
        util::integrate([](double x) { return std::exp(-x); }, 0., 1.);
    assertion(approx_eq(integral, 1. - 1. / M_E),
              std::string{"Numeric integration not correct! Error ~ "} +
                  std::to_string(1. - 1. / M_E - integral));

    double integral2 = util::integrate_coarse(
        [](double x) { return x * x * x + x * x + x - 1; }, 0., 1.);
    assertion(approx_eq(integral2, 1. / 12),
              std::string{"Numeric integration not correct! Error ~ "} +
                  std::to_string(1. / 12 - integral));

    return assertion.status();
}
