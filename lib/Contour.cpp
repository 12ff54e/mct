#include "Contour.h"
#include "util.h"

/**
 * @brief Calculate contour points by finding coordinate between magnetic axis
 * and boundary point that flux(coord) == psi. Thus contour points distributes
 * same as boundary contour points along flux surface.
 *
 */
Contour::Contour(
    double psi,
    const intp::InterpolationFunction<double, 2, MagneticEquilibrium::ORDER>&
        flux,
    const GFileRawData& g_file)
    : flux_(psi), g_file_(g_file) {
    pts_.reserve(g_file.boundary.size());
    for (size_t i = 0; i < g_file.boundary.size(); ++i) {
        pts_.emplace_back(util::vec_field_find_root(flux, g_file.magnetic_axis,
                                                    g_file.boundary[i], psi));
    }

#ifdef _DEBUG
    size_t count{};
    for (size_t i = 0; i < g_file.boundary.size(); ++i) {
        auto pt0 = pts_[i];
        auto pt1 = pts_[(i + 1) % pts_.size()];
        if (util::abs((flux(.5 * (pt0 + pt1)) - psi) / psi) > 1.e-3) {
            ++count;
        }
    }
    if (count > 0) {
        std::cout << "Contour \\psi = " << psi << " has " << count << "/"
                  << pts_.size() << " not-so-accurate segments\n";
    }
#endif
}

size_t Contour::size() const noexcept {
    return pts_.size();
}

double Contour::flux() const noexcept {
    return flux_;
}

const Vec<2, double>& Contour::operator[](size_t i) const {
    return pts_[i];
}
