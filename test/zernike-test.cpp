#include <cmath>
#include <iomanip>
#include <iostream>

#include "BSplineInterpolation/src/include/Mesh.hpp"

#define MCT_MAX_ZERNIKE_ORDER 5
#define MCT_ZERNIKE_POLYNOMIAL_INSTANTIATION
#include "../lib/Zernike.h"

using z_r = std::ratio<0>;
using zernike_radial_5_3 =
    Polynomial<z_r, z_r, z_r, std::ratio<-4>, z_r, std::ratio<5>>;
using zernike_radial_5_3_d =
    Polynomial<z_r, z_r, std::ratio<-12>, z_r, std::ratio<25>>;

int main() {
    double diff{};
    double diff_d{};
    for (int i = 0; i < 10; ++i) {
        double r = .1 * i;
        diff += Zernike::radial_at(5, 3, r) - zernike_radial_5_3::eval(r);
        diff_d += zernike_radial_5_3::derivative<1>(r) -
                  zernike_radial_5_3_d::eval(r);
    }

    std::cout
        << "[" << std::boolalpha << (std::fpclassify(diff) == FP_ZERO)
        << "] Radial part of Zernike basic generated successful (sampling "
           "inspection)\n";
    std::cout << "[" << std::boolalpha << (std::fpclassify(diff_d) == FP_ZERO)
              << "] Derivative of Radial part of Zernike basic generated "
                 "successful (sampling "
                 "inspection)\n";
    std::cout << "["
              << (std::fpclassify(zernike_radial_5_3::derivative<6>(.5)) ==
                  FP_ZERO)
              << "] Derivative order larger than polynomial order test\n";

    constexpr std::size_t nr = 128;
    constexpr std::size_t nt = 128;
    intp::Mesh<double, 2> val(nr, nt);
    std::vector<double> rg(nr);

    constexpr double dr = 1. / static_cast<double>(nr);
    constexpr double dt = 2. * M_PI / static_cast<double>(nt);
    for (std::size_t i = 0; i < nr; ++i) {
        rg[i] = static_cast<double>(i + 1) * dr;
        for (std::size_t j = 0; j < nt; ++j) {
            const double theta = dt * static_cast<double>(j);
            const double r2 = rg[i] * rg[i];
            const double r3 = r2 * rg[i];
            // val = 3*Z(0，0) + Z(2,0) - 1.5*Z(2,2) + 2*Z(5,-3)
            val(i, j) = 2. * (1. + r2) - 1.5 * r2 * std::cos(2. * theta) +
                        2. * r3 * (-4. + 5. * r2) * (std::sin(3. * theta));
        }
    }

    Zernike::Series<double> zernike_series(5, nr, nt, val, rg);

    std::cout << "\nCoefficients of f = 3*Z(0，0) + Z(2,0) - 1.5*Z(2,2) + "
                 "2*Z(5,-3)\n";
    int l = 0;
    for (auto c : zernike_series.coefficient()) {
        std::cout << std::setw(9) << std::setprecision(2) << c << ", ";
        const auto n = Zernike::index_n(l);
        const auto m = Zernike::index_m(l);
        if (n == m) { std::cout << '\n'; }
        ++l;
    }
    return 0;
}
