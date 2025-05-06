#ifndef MEQ_CONTOUR_H
#define MEQ_CONTOUR_H

#include <algorithm>
#include <vector>

#include "GFileRawData.h"
#include "MagneticEquilibrium.h"
#include "Vec.h"
#include "util.h"

class Contour {
   public:
    using pt_type = Vec<2, double>;

   private:
    const double flux_;
    std::vector<pt_type> pts_;
    const GFileRawData& g_file_;

   public:
    Contour(double,
            const intp::
                InterpolationFunction<double, 2, MagneticEquilibrium::ORDER>&,
            const GFileRawData&);

    // properties

    size_t size() const noexcept;

    double flux() const noexcept;

    // element access

    const Vec<2, double>& operator[](size_t) const;

    // numerics method

    template <typename Field, typename Measure>
    double definite_integrate_along(Field, Measure);

    template <typename Field>
    intp::InterpolationFunction1D<3, double> indefinite_integrate_along(
        Field,
        bool = true);
};

template <typename Field, typename Measure>
double Contour::definite_integrate_along(Field f, Measure s) {
    // sort contour pts according to measure
    std::sort(
        pts_.begin(), pts_.end(),
        [&](const pt_type& p1, const pt_type& p2) { return s(p1) < s(p2); });

    std::vector<double> abscissa;
    std::vector<double> ordinate;
    abscissa.reserve(size());
    ordinate.reserve(size());
    for (size_t i = 0; i < size(); ++i) {
        auto& pt = operator[](i);
        abscissa.emplace_back(s(pt.x(), pt.y()));
        ordinate.emplace_back(f(pt.x(), pt.y()));
    }

    // last point being identical to the first one required for periodic
    // interpolation

    abscissa.emplace_back(abscissa.front() + 2 * M_PI);
    ordinate.emplace_back(ordinate.front());

    intp::InterpolationFunction<double, 1, 3> f_interp(
        true, std::make_pair(abscissa.begin(), abscissa.end()),
        std::make_pair(ordinate.begin(), ordinate.end()));

    return util::integrate(f_interp, 0., 2 * M_PI);
}

template <typename Field>
intp::InterpolationFunction1D<3, double> Contour::indefinite_integrate_along(
    Field f,
    bool normalization) {
    std::vector<double> ordinate;
    ordinate.reserve(size() + 1);

    for (auto& pt : pts_) { ordinate.emplace_back(f(pt)); }

    // last point being identical to the first one required for periodic
    // interpolation

    ordinate.push_back(ordinate.front());
    std::vector<double> abscissa;
    if (std::fpclassify(g_file_.geometric_poloidal_angles.front()) == FP_ZERO) {
        abscissa.reserve(size() + 1);
    } else {
        abscissa.reserve(size() + 2);
        abscissa.push_back(0);
    }
    abscissa.insert(abscissa.end(), g_file_.geometric_poloidal_angles.begin(),
                    g_file_.geometric_poloidal_angles.end());
    abscissa.push_back(g_file_.geometric_poloidal_angles.front() + 2 * M_PI);

    intp::InterpolationFunction1D<3, double> f_interp(
        std::make_pair(std::next(abscissa.begin()), abscissa.end()),
        std::make_pair(ordinate.begin(), ordinate.end()), true);

    // do the integration segment by segment

    abscissa.back() = 2 * M_PI;
    std::vector<double> integral;
    integral.reserve(abscissa.size());
    integral.push_back(0);
    for (size_t i = 0; i < abscissa.size() - 1; ++i) {
        integral.push_back(
            util::integrate_coarse(f_interp, abscissa[i], abscissa[i + 1]) +
            integral.back());
    }

    // normalization

    if (normalization) {
        const double coef = 2 * M_PI / integral.back();
        for (auto& v : integral) { v *= coef; }
    }

    intp::InterpolationFunction1D<3, double> integral_interp(
        std::make_pair(abscissa.begin(), abscissa.end()),
        std::make_pair(integral.begin(), integral.end()), true);

    return integral_interp;
}

#endif  // MEQ_CONTOUR_H
