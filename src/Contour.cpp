#include "include/Contour.hpp"
#include "include/util.hpp"

/**
 * @brief Calculate contour points by finding coordinate between magnetic axis
 * and boundary point that flux(coord) == psi. Thus contour points distributes
 * same as boundary contour points along flux surface.
 *
 * @param psi
 * @param flux
 * @param g_file
 */
Contour::Contour(double psi,
                 const intp::InterpolationFunction<double, 2>& flux,
                 const GFileRawData& g_file)
    : g_file(g_file) {
    pts.reserve(g_file.boundary.size());
    for (size_t i = 0; i < g_file.boundary.size(); ++i) {
        pts.emplace_back(util::vec_field_find_root(flux, g_file.magnetic_axis,
                                                   g_file.boundary[i], psi));
    }
}

size_t Contour::size() const noexcept {
    return pts.size();
}

const Vec<2, double>& Contour::operator[](size_t i) const {
    return pts[i];
}