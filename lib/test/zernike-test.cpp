#include <cmath>
#include <iomanip>
#include <iostream>

#include "BSplineInterpolation/src/include/Mesh.hpp"

#define ZQ_TIMER_IMPLEMENTATION
#include "Timer.h"

#define MEQ_MAX_ZERNIKE_POLAR_ORDER 6
#define MEQ_ZERNIKE_POLYNOMIAL_INSTANTIATION
#include "Zernike.h"

using zernike_radial_5_3 = Polynomial<0, 0, 0, -4, 0, 5>;
using zernike_radial_5_3_d = Polynomial<0, 0, -12, 0, 25>;

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

    constexpr std::size_t nr = 256;
    constexpr std::size_t nt = 256;
    intp::Mesh<double, 2> val(nr, nt);
    std::vector<double> rg(nr);

    auto& timer = Timer::get_timer();
    timer.start("Grid");

    constexpr double dr = 1. / static_cast<double>(nr);
    constexpr double dt = 2. * M_PI / static_cast<double>(nt);
    for (std::size_t i = 0; i < nr; ++i) {
        rg[i] = static_cast<double>(i + 1) * dr;
        const double r2 = rg[i] * rg[i];
        const double r3 = r2 * rg[i];
        for (std::size_t j = 0; j < nt; ++j) {
            const double theta = dt * static_cast<double>(j);
            // val = 3*Z(0，0) + Z(2,0) - 1.5*Z(2,2) + 2*Z(5,-3)
            val(i, j) = 3. + (2. * r2 - 1.) - 1.5 * r2 * std::cos(2. * theta) +
                        2. * r3 * (-4. + 5. * r2) * (std::sin(3. * theta));
        }
    }

    timer.pause_last_and_start_next("Series");

    constexpr int polar_order = 4;
    Zernike::Series<double> zernike_series(polar_order, nr, nt, val, rg);

    timer.pause();
    std::cout << "\nCoefficients of f = 3*Z(0，0) + Z(2,0) - 1.5*Z(2,2) + "
                 "2*Z(5,-3)\n";
    int l = 0;

    const int default_precision = static_cast<int>(std::cout.precision(2));
    for (auto c : zernike_series.coefficient()) {
        std::cout << std::setw(9) << c << ", ";
        const auto [n, m] = Zernike::basic_index_nm(l, zernike_series.order);
        if (n == m || n + m == 2 * static_cast<int>(zernike_series.order)) {
            std::cout << '\n';
        }
        ++l;
    }
    std::cout << std::setprecision(default_precision);

    timer.start("Evaluate");

    diff = 0;
    for (std::size_t i = 0; i < nr; ++i) {
        const auto r = rg[i];
        for (std::size_t j = 0; j < nt; ++j) {
            const double theta = dt * static_cast<double>(j);
            diff += std::pow(val(i, j) - zernike_series(r, theta), 2);
        }
    }

    timer.pause();
    std::cout << "Grid Size: " << nr << " * " << nt << '\n'
              << "L2 Difference on Grid: " << std::sqrt(diff / nr * nt) << '\n';

    timer.print();

    return 0;
}
