#include "include/Contour.hpp"
#include "include/util.hpp"

/**
 * @brief Calculate contour points by finding coordinate between magnetic axis
 * and boundary point that flux(coord) == psi. Thus contour points distributes
 * same as boundary contour points along flux surface.
 *
 * @param psi
 * @param flux
 * @param gfile
 */
Contour::Contour(double psi,
                 const intp::InterpolationFunction<double, 2>& flux,
                 const GFileRawData& gfile)
    : gfile(gfile) {
    pts.reserve(gfile.boundary.size());
    for (size_t i = 0; i < gfile.boundary.size(); ++i) {
        pts.emplace_back(util::vec_field_find_root(flux, gfile.magnetic_axis,
                                                   gfile.boundary[i], psi));
    }
}

size_t Contour::size() const noexcept {
    return pts.size();
}

const Vec<2, double>& Contour::operator[](size_t i) const {
    return pts[i];
}