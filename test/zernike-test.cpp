#include <cmath>
#include <iomanip>
#include <iostream>

#define MCT_MAX_ZERNIKE_ORDER 5
#include "../lib/Zernike.h"

using z_r = std::ratio<0>;
using zernike_radial_5_3 =
    Polynomial<z_r, z_r, z_r, std::ratio<-4>, z_r, std::ratio<5>>;

int main() {
    double diff = 0;
    for (int i = 0; i < 10; ++i) {
        double r = .1 * i;
        diff += zernike_radial_at(5, 3, r) - zernike_radial_5_3::eval(r);
    }

    std::cout << "Radial part of Zernike basic generated successful (sampling "
                 "inspection): "
              << std::boolalpha << (std::fpclassify(diff) == FP_ZERO) << '\n';
    return 0;
}
