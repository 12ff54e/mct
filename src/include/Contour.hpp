#pragma once

#include <algorithm>
#include <vector>

#include "lib/BSplineInterpolation/intp"

#include "GFileRawData.hpp"
#include "Vec.hpp"
#include "util.hpp"

class Contour {
   public:
    using pt_type = Vec<2, double>;

   private:
    std::vector<pt_type> pts;
    const GFileRawData& gfile;

   public:
    Contour(double,
            const intp::InterpolationFunction<double, 2>&,
            const GFileRawData&);

    // properties

    size_t size() const noexcept;

    // elment access

    const Vec<2, double>& operator[](size_t) const;

    // numerics method

    template <typename Field, typename Measure>
    double definite_integrate_along(Field, Measure);

    template <typename Field>
    intp::InterpolationFunction1D<double> indefinite_integrate_along(
        Field,
        bool = true);
};

template <typename Field, typename Measure>
double Contour::definite_integrate_along(Field f, Measure s) {
    // sort contour pts according to measure
    std::sort(
        pts.begin(), pts.end(),
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

    intp::InterpolationFunction<double, 1> f_interp(
        3, true, std::make_pair(abscissa.begin(), abscissa.end()),
        std::make_pair(ordinate.begin(), ordinate.end()));

    return util::integrate(f_interp, 0., 2 * M_PI);
}

template <typename Field>
intp::InterpolationFunction1D<double> Contour::indefinite_integrate_along(
    Field f,
    bool normalization) {
    std::vector<double> ordinate;
    ordinate.reserve(size() + 1);

    for (auto& pt : pts) { ordinate.emplace_back(f(pt)); }

    // last point being identical to the first one required for periodic
    // interpolation

    ordinate.push_back(ordinate.front());
    std::vector<double> abscissa;
    if (gfile.geometric_poloidal_angles.front() == 0) {
        abscissa.reserve(size() + 1);
    } else {
        abscissa.reserve(size() + 2);
        abscissa.push_back(0);
    }
    abscissa.insert(abscissa.end(), gfile.geometric_poloidal_angles.begin(),
                    gfile.geometric_poloidal_angles.end());
    abscissa.push_back(gfile.geometric_poloidal_angles.front() + 2 * M_PI);

    intp::InterpolationFunction1D<double> f_interp(
        std::make_pair(std::next(abscissa.begin()), abscissa.end()),
        std::make_pair(ordinate.begin(), ordinate.end()), 3, true);

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

    intp::InterpolationFunction1D<double> integral_interp(
        std::make_pair(abscissa.begin(), abscissa.end()),
        std::make_pair(integral.begin(), integral.end()), 3, true);

    return integral_interp;
}
