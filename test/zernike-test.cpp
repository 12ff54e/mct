#include <cmath>
#include <iomanip>
#include <iostream>

#define MCT_MAX_ZERNIKE_ORDER 5
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
    return 0;
}
