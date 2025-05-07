#include <iomanip>
#include <limits>
#include <stdexcept>  // runtime_error

#include "Contour.h"
#include "MagneticEquilibrium.h"
#include "Timer.h"

#ifdef MEQ_ZERNIKE_SERIES_
#define MEQ_ZERNIKE_POLYNOMIAL_INSTANTIATION
#include "Zernike.h"
#endif

MagneticEquilibrium::MagneticEquilibriumRaw_
MagneticEquilibrium::generate_boozer_coordinate_(
    const GFileRawData& g_file_data,
    std::size_t radial_sample,
    double psi_ratio) {
    if (radial_sample % 2 != 0) {
        throw std::runtime_error(
            "[MagneticEquilibrium] Radial sample point must be even.");
    }

    auto& timer = Timer::get_timer();
    timer.start("Create Boozer grid");
    intp::InterpolationFunction<double, 2u, ORDER> flux_function(
        g_file_data.flux,
        std::make_pair(g_file_data.r_left,
                       g_file_data.r_left + g_file_data.dim.x()),
        std::make_pair(g_file_data.z_mid - .5 * g_file_data.dim.y(),
                       g_file_data.z_mid + .5 * g_file_data.dim.y()));

    // check poloidal flux value at magnetic axis
    const auto psi_ma_intp = flux_function(g_file_data.magnetic_axis);
    if (std::abs((psi_ma_intp - g_file_data.flux_magnetic_axis) /
                 (g_file_data.flux_LCFS - g_file_data.flux_magnetic_axis)) >
        1.e-4) {
        std::cout << "The poloidal flux of magnetic axis given in gfile "
                     "deviates from interpolated value.\n"
                  << "  \\psi_p in gfile: " << g_file_data.flux_magnetic_axis
                  << "\n  \\psi_p from interpolation: " << psi_ma_intp << '\n';
    }

    double psi_boundary_min = 10. * std::pow(g_file_data.r_center, 2) *
                              std::abs(g_file_data.b_center);
    for (const auto& pt : g_file_data.boundary) {
        psi_boundary_min = std::min(psi_boundary_min, flux_function(pt));
    }
    std::cout << "The poloidal flux of last closed flux surface is "
              << g_file_data.flux_LCFS << '\n'
              << "Minimum of interpolated value at boundary points is "
              << psi_boundary_min << '\n';

    const double psi_bd = g_file_data.flux_LCFS - psi_ma_intp;
    double psi_wall = psi_ratio * psi_bd;
    if (psi_wall > psi_boundary_min - psi_ma_intp) {
        psi_wall = psi_boundary_min - psi_ma_intp;
        std::cout << "Interpolated flux value at boundary is too small, so "
                     "psi_wall is set to this value.\n";
    }
    psi_delta_ = psi_wall / static_cast<double>(lsp - 1);

    // contours are from \\Delta\\psi to LCFS
    std::vector<Contour> contours;
    contours.reserve(radial_sample);
    for (std::size_t i = 0; i < radial_sample; ++i) {
        contours.emplace_back(
            util::lerp(psi_delta_, psi_wall,
                       static_cast<double>(i) /
                           static_cast<double>(radial_sample - 1)) +
                psi_ma_intp,
            flux_function, g_file_data);
    }

    constexpr double magnetic_constant = 4.e-7 * M_PI;

    // safety factor, on shifted psi
    intp::InterpolationFunction1D<ORDER> safety_factor_intp{
        std::make_pair(0., psi_bd),
        intp::util::get_range(g_file_data.safety_factor)};

    // poloidal current, on shifted psi
    intp::InterpolationFunction1D<ORDER> poloidal_current_intp{
        std::make_pair(0., psi_bd), intp::util::get_range(g_file_data.f_pol)};

    // pressure, on shifted psi
    intp::InterpolationFunction1D<ORDER> pressure_intp{
        std::make_pair(0., psi_bd),
        intp::util::get_range(g_file_data.pressure)};

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
    intp::InterpolationFunctionTemplate1D<ORDER> poloidal_template{
        intp::util::get_range(poloidal_angles), poloidal_angles.size(), true};

    if (std::fpclassify(poloidal_angles.front()) != FP_ZERO) {
        poloidal_angles.insert(poloidal_angles.begin(), 0);
    }

    poloidal_angles.back() = PI2;
    // \\theta range: 0, ..., 2\\pi
    intp::InterpolationFunctionTemplate1D<ORDER> poloidal_template_full{
        intp::util::get_range(poloidal_angles), poloidal_angles.size(), false};

    // output data

    intp::Mesh<double, 2> magnetic_boozer(radial_sample, lst + 1);
    intp::Mesh<double, 2> r_boozer(radial_sample, lst + 1);
    intp::Mesh<double, 2> z_boozer(radial_sample, lst + 1);
    intp::Mesh<double, 2> jacobian_boozer(radial_sample, lst + 1);

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
        const double psi = contours[ri].flux() - psi_ma_intp;
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
        const auto b2j_flux_avg = b2j_int.back() / PI2;
        tor_current_n.push_back(bp2j_int.back() / (PI2 * current_unit));
        // normalization
        for (auto& v : b2j_int) { v /= b2j_flux_avg; }
        auto boozer_geo_intp =
            poloidal_template_full.interpolate(intp::util::get_range(b2j_int));

        // calculate necessary values on a even-spaced boozer grid
        for (size_t i = 0; i <= lst; ++i) {
            double theta_boozer =
                (static_cast<double>(i % lst) + .5) * theta_delta_;
            double theta_geo = util::find_root(
                [&](double t) { return boozer_geo_intp(t) - theta_boozer; }, 0.,
                PI2);
            double r_grid = r_geo_intp(theta_geo);
            double z_grid = z_geo_intp(theta_geo);

            // be careful of normalization
            const auto b = b_field({r_grid, z_grid}, psi);
            magnetic_boozer(ri, i) = b / magnetic_field_unit;
            r_boozer(ri, i) = r_grid / length_unit;
            // z value is shifted such that magnetic axis has z = 0
            z_boozer(ri, i) =
                (z_grid - g_file_data.magnetic_axis.y()) / length_unit;
            jacobian_boozer(ri, i) =
                b2j_flux_avg / (b * b) * magnetic_field_unit / length_unit;
        }

        safety_factor.push_back(safety_factor_intp(psi));
        pol_current_n.push_back(poloidal_current_intp(psi) / current_unit);
        pressure_n.push_back(pressure_intp(psi) / pressure_unit);
        tor_flux_n.push_back(
            (ri == 0 ? 0. : tor_flux_n.back()) +
            util::integrate_coarse(
                safety_factor_intp,
                ri == 0 ? 0. : (contours[ri - 1].flux() - psi_ma_intp), psi) /
                flux_unit);
        // r_minor defined as distance from magnetic axis at weak field side
        // this value is always normalized to R0
        r_minor_n.push_back(r_geo_intp(0.) / R0 - 1.);
    }

    timer.pause();

    // psi_delta_ is normalized after flux surface is fully constructed, and
    // should never be changed hereafter
    psi_delta_ /= flux_unit;

    const double q0 = safety_factor_intp(0);
    const double b0n = B0 / magnetic_field_unit;
    const double g0n = poloidal_current_intp(0) / current_unit;
    const double p0n = pressure_intp(0) / pressure_unit;
    return MagneticEquilibriumRaw_{
        {std::move(magnetic_boozer), std::move(r_boozer), std::move(z_boozer),
         std::move(jacobian_boozer)},
        {std::move(safety_factor), std::move(pol_current_n),
         std::move(tor_current_n), std::move(pressure_n), std::move(r_minor_n),
         std::move(tor_flux_n)},
        {b0n, R0 / length_unit, 0., q0 * g0n / (b0n * b0n)},
        {q0, g0n, 0., p0n, 0., 0.},
        flux_unit};
}

double MagneticEquilibrium::psi_delta() const {
    return psi_delta_;
}
double MagneticEquilibrium::theta_delta() const {
    return theta_delta_;
}

const std::array<double, MagneticEquilibrium::FIELD_NUM_2D>&
MagneticEquilibrium::axis_value_2d() const {
    return spdata_raw_.axis_value_2d;
}
const std::array<double, MagneticEquilibrium::FIELD_NUM_1D>&
MagneticEquilibrium::axis_value_1d() const {
    return spdata_raw_.axis_value_1d;
}
const std::vector<double>& MagneticEquilibrium::psi_for_output() const {
    return intp_data().psi_sample_for_output;
}

const MagneticEquilibrium::MagneticEquilibriumIntp_&
MagneticEquilibrium::intp_data() const {
    return spdata_intp_;
}

#ifdef MEQ_ZERNIKE_SERIES_
Zernike::Series<double>
#else
intp::InterpolationFunction<double, 2, MagneticEquilibrium::ORDER_OUT>
#endif
MagneticEquilibrium::create_2d_spline_(
    const intp::Mesh<double, 2>& data,
    const std::vector<double>& psi_sample) const {
#ifdef MEQ_ZERNIKE_SERIES_
    static_cast<void>(psi_sample);
    // The Zernike series is actually representing f(r, theta+delta/2)

    std::vector<double> r(spdata_raw_.data_1d[5]);
    const auto psi_w = r[r.size() - 1];
    for (auto& v : r) { v = std::sqrt(v / psi_w); }

    const auto zernike_order = static_cast<int>(
        lst / 5 > MEQ_MAX_ZERNIKE_POLAR_ORDER ? MEQ_MAX_ZERNIKE_POLAR_ORDER
                                              : lst / 5);
    return {zernike_order, r.size(), lst, data, r};
#else
    // interpolate the even-spaced data
    intp::InterpolationFunction<double, 2, ORDER_OUT> data_intp(
        {false, true}, data,
        std::make_pair(psi_delta(), psi_delta() * static_cast<double>(lsp - 1)),
        std::make_pair(.5 * theta_delta_, 2. * M_PI + .5 * theta_delta_));

    if (psi_sample.empty()) { return data_intp; }

    // resample on the interpolated function
    intp::Mesh<double, 2> data_resampled(lsp, lst + 1);
    for (std::size_t i = 0; i < lsp; ++i) {
        for (std::size_t j = 0; j <= lst; ++j) {
            data_resampled(i, j) = data_intp(
                psi_sample[i], (static_cast<double>(j) + .5) * theta_delta_);
        }
    }
    // construct the interpolation function for output
    return intp::InterpolationFunction<double, 2, ORDER_OUT>{
        {false, true},
        std::move(data_resampled),
        intp::util::get_range(psi_sample),
        std::make_pair(.5 * theta_delta_, 2. * M_PI + .5 * theta_delta_)};
#endif
}

intp::InterpolationFunction1D<MagneticEquilibrium::ORDER_OUT, double>
MagneticEquilibrium::create_1d_spline_(
    const std::vector<double>& data,
    const std::vector<double>& psi_sample) const {
    // interpolate the even-spaced data
    intp::InterpolationFunction1D<ORDER_OUT, double> data_intp(
        std::make_pair(psi_delta(), psi_delta() * static_cast<double>(lsp - 1)),
        intp::util::get_range(data), false);

    if (psi_sample.empty()) { return data_intp; }

    // resample on the interpolated function
    std::vector<double> data_resampled;
    data_resampled.reserve(psi_sample.size());
    for (const auto& psi : psi_sample) {
        data_resampled.push_back(data_intp(psi));
    }

    return intp::InterpolationFunction1D<ORDER_OUT, double>{
        intp::util::get_range(psi_sample),
        intp::util::get_range(data_resampled)};
}
