#include <iomanip>
#include <limits>

#include "include/Contour.hpp"
#include "include/Spdata.hpp"

Spdata::Spdata(const GFileRawData& g_file_data,
               std::size_t radial_grid_num,
               std::size_t poloidal_grid_num,
               bool use_si,
               std::size_t radial_sample,
               double psi_ratio)
    : use_si_(use_si),
      lsp_(radial_grid_num),
      lst_(poloidal_grid_num),
      psi_delta_(psi_ratio *
                 (g_file_data.flux_LCFS - g_file_data.flux_magnetic_axis) /
                 static_cast<double>(lsp_ - 1)),
      theta_delta_(2. * M_PI / static_cast<double>(lst_)),
      spdata_raw_{
          generate_boozer_coordinate_(g_file_data, radial_sample, psi_ratio)},
      spdata_intp_{spdata_raw_, *this, std::make_index_sequence<FIELD_NUM_2D>{},
                   std::make_index_sequence<FIELD_NUM_1D>{}} {}

std::ostream& operator<<(std::ostream& os, const Spdata& spdata) {
    // mesh grid info
    os << std::setw(4) << spdata.lsp_ << std::setw(4) << spdata.lst_
       << std::setw(4) << 4 << std::setw(4) << 8 << '\n';
    os << std::scientific << std::uppercase << std::setprecision(10);
    os << std::setw(18) << spdata.spdata_intp_.psi_sample_for_output.back()
       << std::setw(18) << spdata.spdata_intp_.psi_sample_for_output.back()
       << '\n';

    const auto psi_delta = spdata.psi_delta_ / spdata.spdata_raw_.flux_unit;
    const auto theta_delta = spdata.theta_delta_;
    const auto lst = spdata.lst_;

    auto write_1d_coef = [&](const auto& f_1d, std::size_t idx,
                             double value_on_axis, bool singular) {
        double c0, c1, c2;
        if (idx == 0) {
            const double v1 = f_1d(psi_delta) - value_on_axis;
            const double v2 = f_1d.derivative(std::make_pair(psi_delta, 1)) *
                              (singular ? 2. * std::sqrt(psi_delta) : 1.);
            const double d = singular ? std::sqrt(psi_delta) : psi_delta;
            c0 = value_on_axis;
            c1 = 2. * v1 / d - v2;
            c2 = (-v1 + v2 * d) / (d * d);
        } else {
            const double psi = static_cast<double>(idx) * psi_delta;
            c0 = f_1d(psi);
            c1 = f_1d.derivative(std::make_pair(psi, 1));
            c2 = .5 * f_1d.derivative(std::make_pair(psi + .5 * psi_delta, 2));
        }
        os << std::setw(18) << c0 << std::setw(18) << c1 << std::setw(18) << c2
           << '\n';
    };

    auto write_2d_coef = [&](auto& f_2d, std::size_t idx,
                             double value_on_axis) {
        for (size_t i = 0; i < 9; ++i) {
            size_t psi_order = i % 3;
            size_t theta_order = i / 3;
            for (size_t j = 0; j <= lst; ++j) {
                double coef =
                    (psi_order == 2 ? .5 : 1.) * (theta_order == 2 ? .5 : 1.);
                const double theta =
                    (static_cast<double>(j) + (theta_order == 2 ? .5 : 0.)) *
                    theta_delta;
                if (idx == 0) {
                    const double v1 =
                        theta_order == 0
                            ? f_2d(psi_delta, theta) - value_on_axis
                            : f_2d.derivative({psi_delta, theta},
                                              {0, theta_order});
                    const double v2 =
                        f_2d.derivative({psi_delta, theta}, {1, theta_order});

                    coef = psi_order == 0
                               ? theta_order == 0 ? value_on_axis : 0.
                               : 2. * coef *
                                     (psi_order == 1
                                          ? (v1 - v2 * psi_delta) /
                                                std::sqrt(psi_delta)
                                          : -(v1 - 2. * v2 * psi_delta) /
                                                psi_delta);
                } else {
                    const double psi = (static_cast<double>(idx) +
                                        (psi_order == 2 ? .5 : 0.)) *
                                       psi_delta;
                    coef = (i == 0 ? f_2d(psi, theta)
                                   : coef * f_2d.derivative(
                                                {psi, theta},
                                                {psi_order, theta_order}));
                }

                os << std::setw(18) << coef;
                if (j % 4 == 3) { os << '\n'; }
            }
            if ((lst + 1) % 4 != 0) { os << '\n'; }
        }
    };

    for (std::size_t idx = 0; idx < spdata.lsp_; ++idx) {
        for (std::size_t i_2d = 0; i_2d < Spdata::FIELD_NUM_2D; ++i_2d) {
            write_2d_coef(spdata.spdata_intp_.intp_2d[i_2d], idx,
                          spdata.spdata_raw_.axis_value_2d[i_2d]);
        }
        for (std::size_t i_1d = 0; i_1d < Spdata::FIELD_NUM_1D; ++i_1d) {
            write_1d_coef(spdata.spdata_intp_.intp_1d[i_1d], idx,
                          spdata.spdata_raw_.axis_value_1d[i_1d], i_1d == 4);
        }
    }

    // TODO: ripple related
    os << std::setw(4) << 0 << std::setw(4) << 0 << '\n';
    os << std::setw(18) << spdata.spdata_raw_.axis_value_2d[1] << std::setw(18)
       << 0. << std::setw(18) << 0. << '\n';
    os << std::setw(18) << 0. << std::setw(18) << 0. << '\n';

    return os;
}

Spdata::SpdataRaw_ Spdata::generate_boozer_coordinate_(
    const GFileRawData& g_file_data,
    std::size_t radial_sample,
    double psi_ratio) {
    intp::InterpolationFunction<double, 2u> flux_function(
        ORDER_, g_file_data.flux,
        std::make_pair(g_file_data.r_center - .5 * g_file_data.dim.x(),
                       g_file_data.r_center + .5 * g_file_data.dim.x()),
        std::make_pair(g_file_data.z_mid - .5 * g_file_data.dim.y(),
                       g_file_data.z_mid + .5 * g_file_data.dim.y()));

    const double psi_bd =
        g_file_data.flux_LCFS - g_file_data.flux_magnetic_axis;
    const double psi_wall = psi_ratio * psi_bd;

    // contours are from \\Delta\\psi to LCFS
    std::vector<Contour> contours;
    contours.reserve(radial_sample);
    for (std::size_t i = 0; i < radial_sample; ++i) {
        contours.emplace_back(
            util::lerp(psi_delta_, psi_wall,
                       static_cast<double>(i) /
                           static_cast<double>(radial_sample - 1)) +
                g_file_data.flux_magnetic_axis,
            flux_function, g_file_data);
    }

    constexpr double magnetic_constant = 4.e-7 * M_PI;

    // safety factor, on shifted psi
    intp::InterpolationFunction1D<> safety_factor_intp{
        std::make_pair(0., psi_bd),
        intp::util::get_range(g_file_data.safety_factor), ORDER_};

    // poloidal current, on shifted psi
    intp::InterpolationFunction1D<> poloidal_current_intp{
        std::make_pair(0., psi_bd), intp::util::get_range(g_file_data.f_pol),
        ORDER_};

    // pressure, on shifted psi
    intp::InterpolationFunction1D<> pressure_intp{
        std::make_pair(0., psi_bd), intp::util::get_range(g_file_data.pressure),
        ORDER_};

    // this following function accepts shifted psi (0 at m.a.)

    auto b2j_field = [&](Vec<2, double> pt, double psi) {
        double dp_dr = flux_function.derivative(pt, {1, 0});
        double dp_dz = flux_function.derivative(pt, {0, 1});
        double r = pt.x();
        pt -= g_file_data.magnetic_axis;
        double r2 = pt.L2_norm_square_();
        double f = poloidal_current_intp(psi);

        return (f * f + dp_dr * dp_dr + dp_dz * dp_dz) * r2 /
               (r * (dp_dr * pt.x() + dp_dz * pt.y()));
    };

    auto b_field = [&](Vec<2, double> pt, double psi) {
        double dp_dr = flux_function.derivative(pt, {1, 0});
        double dp_dz = flux_function.derivative(pt, {0, 1});
        double f = poloidal_current_intp(psi);

        return std::sqrt(f * f + dp_dr * dp_dr + dp_dz * dp_dz) / pt.x();
    };

    auto j_field = [&](Vec<2, double> pt, double) {
        double dp_dr = flux_function.derivative(pt, {1, 0});
        double dp_dz = flux_function.derivative(pt, {0, 1});
        double r = pt.x();
        pt -= g_file_data.magnetic_axis;
        double r2 = pt.L2_norm_square_();

        return r * r2 / (dp_dr * pt.x() + dp_dz * pt.y());
    };

    auto bp2j_field = [&](Vec<2, double> pt, double) {
        double dp_dr = flux_function.derivative(pt, {1, 0});
        double dp_dz = flux_function.derivative(pt, {0, 1});
        double r = pt.x();
        pt -= g_file_data.magnetic_axis;
        double r2 = pt.L2_norm_square_();

        return (dp_dr * dp_dr + dp_dz * dp_dz) * r2 /
               ((dp_dr * pt.x() + dp_dz * pt.y()) * r);
    };

    constexpr double PI2 = 2 * M_PI;

    std::vector<double> poloidal_angles{g_file_data.geometric_poloidal_angles};
    poloidal_angles.push_back(poloidal_angles.front() + PI2);
    // \\theta range: \\theta_0, ..., \\theta_0 + 2\\pi
    intp::InterpolationFunctionTemplate1D<> poloidal_template{
        intp::util::get_range(poloidal_angles), poloidal_angles.size(), 5,
        true};

    if (std::fpclassify(poloidal_angles.front()) != FP_ZERO) {
        poloidal_angles.insert(poloidal_angles.begin(), 0);
    }

    poloidal_angles.back() = PI2;
    // \\theta range: 0, ..., 2\\pi
    intp::InterpolationFunctionTemplate1D<> poloidal_template_full{
        intp::util::get_range(poloidal_angles), poloidal_angles.size(), 5,
        false};

    // output data

    intp::Mesh<double, 2> magnetic_boozer(radial_sample, lst_ + 1);
    intp::Mesh<double, 2> r_boozer(radial_sample, lst_ + 1);
    intp::Mesh<double, 2> z_boozer(radial_sample, lst_ + 1);
    intp::Mesh<double, 2> jacobian_boozer(radial_sample, lst_ + 1);

    std::vector<double> safety_factor, pol_current_n, tor_current_n, pressure_n,
        r_minor_n, tor_flux_n;

    safety_factor.reserve(radial_sample);
    pol_current_n.reserve(radial_sample);
    tor_current_n.reserve(radial_sample);
    pressure_n.reserve(radial_sample);
    r_minor_n.reserve(radial_sample);
    tor_flux_n.reserve(radial_sample);

    const double B0 = b_field(g_file_data.magnetic_axis, 0.);
    const double R0 = g_file_data.magnetic_axis.x();

    // This two basic unit determines the output spdata unit,
    // setting them to 1 means SI unit.
    const double length_unit = use_si_ ? 1. : R0;
    const double magnetic_field_unit = use_si_ ? 1. : B0;

    const double current_unit = length_unit * magnetic_field_unit;
    const double pressure_unit =
        magnetic_field_unit * magnetic_field_unit / magnetic_constant;
    const double flux_unit = length_unit * length_unit * magnetic_field_unit;

    // construct boozer coordinate, integrate B^2 * J along each contour

#define boozer_list() \
    X(r);             \
    X(z);             \
    X(b2j);           \
    X(bp2j)

    for (std::size_t ri = 0; ri < contours.size(); ++ri) {
        const double psi = contours[ri].flux() - g_file_data.flux_magnetic_axis;
        const std::size_t poloidal_size = contours[ri].size() + 1;
#define X(name)                     \
    std::vector<double> name##_geo; \
    name##_geo.reserve(poloidal_size)
        boozer_list();
#undef X

        // quantities on geometric grid
        for (size_t i = 0; i < poloidal_size; ++i) {
            const auto& pt = contours[ri][i % (poloidal_size - 1)];
            r_geo.push_back(pt.x());
            z_geo.push_back(pt.y());
            b2j_geo.push_back(b2j_field(pt, psi));
            bp2j_geo.push_back(bp2j_field(pt, psi));
        }
        // interpolation on geometric grid
#define X(name)            \
    auto name##_geo_intp = \
        poloidal_template.interpolate(intp::util::get_range(name##_geo))
        boozer_list();
#undef X

        // integrate

        std::vector<double> b2j_int;
        b2j_int.reserve(poloidal_angles.size());
        b2j_int.push_back(0);

        std::vector<double> bp2j_int;
        bp2j_int.reserve(poloidal_angles.size());
        bp2j_int.push_back(0);

        // Poloidal grid begins from \\theta = 0 and ends at \\theta = 2\\pi
        for (size_t i = 1; i < poloidal_angles.size(); ++i) {
            b2j_int.push_back(b2j_int.back() +
                              util::integrate_coarse(b2j_geo_intp,
                                                     poloidal_angles[i - 1],
                                                     poloidal_angles[i]));
            bp2j_int.push_back(bp2j_int.back() +
                               util::integrate_coarse(bp2j_geo_intp,
                                                      poloidal_angles[i - 1],
                                                      poloidal_angles[i]));
        }
        auto coef = b2j_int.back() / PI2;
        tor_current_n.push_back(bp2j_int.back() / (PI2 * current_unit));
        // normalization
        for (auto& v : b2j_int) { v /= coef; }
        auto boozer_geo_intp =
            poloidal_template_full.interpolate(intp::util::get_range(b2j_int));

        // calculate necessary values on a even-spaced boozer grid
        for (size_t i = 0; i <= lst_; ++i) {
            double theta_boozer =
                (static_cast<double>(i % lst_) + .5) * theta_delta_;
            double theta_geo = util::find_root(
                [&](double t) { return boozer_geo_intp(t) - theta_boozer; }, 0.,
                PI2);
            double r_grid = r_geo_intp(theta_geo);
            double z_grid = z_geo_intp(theta_geo);

            // be careful of normalization

            magnetic_boozer(ri, i) =
                b_field({r_grid, z_grid}, psi) / magnetic_field_unit;
            r_boozer(ri, i) = r_grid / length_unit;
            // z value is shifted such that magnetic axis has z = 0
            z_boozer(ri, i) =
                (z_grid - g_file_data.magnetic_axis.y()) / length_unit;
            jacobian_boozer(ri, i) = j_field({r_grid, z_grid}, psi) *
                                     magnetic_field_unit / length_unit;
        }

        safety_factor.push_back(safety_factor_intp(psi));
        pol_current_n.push_back(poloidal_current_intp(psi) / current_unit);
        pressure_n.push_back(pressure_intp(psi) / pressure_unit);
        tor_flux_n.push_back(
            (ri == 0 ? 0. : tor_flux_n.back()) +
            util::integrate_coarse(safety_factor_intp,
                                   ri == 0 ? 0.
                                           : (contours[ri - 1].flux() -
                                              g_file_data.flux_magnetic_axis),
                                   psi) /
                flux_unit);
        // r_minor defined as distance from magnetic axis at weak field side
        // this value is always normalized to R0
        r_minor_n.push_back(r_geo_intp(0.) / R0 - 1.);
    }
    const double q0 = safety_factor_intp(0);
    const double b0n = B0 / magnetic_field_unit;
    const double g0n = poloidal_current_intp(0) / current_unit;
    const double p0n = pressure_intp(0) / pressure_unit;
    return SpdataRaw_{{std::move(magnetic_boozer), std::move(r_boozer),
                       std::move(z_boozer), std::move(jacobian_boozer)},
                      {std::move(safety_factor), std::move(pol_current_n),
                       std::move(tor_current_n), std::move(pressure_n),
                       std::move(r_minor_n), std::move(tor_flux_n)},
                      {b0n, R0 / length_unit, 0., q0 * g0n / (b0n * b0n)},
                      {q0, g0n, 0., p0n, 0., 0.},
                      flux_unit};
}

std::vector<double> Spdata::generate_psi_sample_for_output_(double unit) const {
    std::vector<double> psi(lsp_);
    psi[0] = psi_delta_ / unit;
    for (std::size_t i = 1; i < lsp_ - 1; ++i) {
        psi[i] = psi_delta_ * (static_cast<double>(i) + .5) / unit;
    }
    psi[lsp_ - 1] = static_cast<double>(lsp_ - 1) * psi_delta_ / unit;

    return psi;
}

intp::InterpolationFunction<double, 2> Spdata::create_2d_spline_(
    const intp::Mesh<double, 2>& data,
    const std::vector<double>& psi_sample) const {
    // interpolate the even-spaced data
    intp::InterpolationFunction<double, 2> data_intp(
        ORDER_, {false, true}, data,
        std::make_pair(psi_sample.front(), psi_sample.back()),
        std::make_pair(.5 * theta_delta_, 2. * M_PI + .5 * theta_delta_));

    // resample on the interpolated function
    intp::Mesh<double, 2> data_resampled(lsp_, lst_ + 1);
    for (std::size_t i = 0; i < lsp_; ++i) {
        for (std::size_t j = 0; j <= lst_; ++j) {
            data_resampled(i, j) = data_intp(
                psi_sample[i], (static_cast<double>(j) + .5) * theta_delta_);
        }
    }
    // construct the interpolation function for output
    return intp::InterpolationFunction<double, 2>{
        ORDER_OUT_,
        {false, true},
        std::move(data_resampled),
        intp::util::get_range(psi_sample),
        std::make_pair(.5 * theta_delta_, 2. * M_PI + .5 * theta_delta_)};
}

intp::InterpolationFunction1D<double> Spdata::create_1d_spline_(
    const std::vector<double>& data,
    const std::vector<double>& psi_sample) const {
    // interpolate the even-spaced data
    intp::InterpolationFunction1D<double> data_intp(
        std::make_pair(psi_sample.front(), psi_sample.back()),
        intp::util::get_range(data), ORDER_, false);

    // resample on the interpolated function
    std::vector<double> data_resampled;
    data_resampled.reserve(psi_sample.size());
    for (const auto& psi : psi_sample) {
        data_resampled.push_back(data_intp(psi));
    }

    return intp::InterpolationFunction1D<double>{
        intp::util::get_range(psi_sample),
        intp::util::get_range(data_resampled), ORDER_OUT_};
}
