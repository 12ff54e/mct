#include <cerrno>
#include <chrono>
#include <fstream>
#include <iomanip>

#include "./include/Contour.hpp"
#include "./include/GFileRawData.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(int argc, char** argv) {
    using namespace std::chrono;

    auto t_start = high_resolution_clock::now();

    if (argc == 1) {
        std::cout << "Please provide gfile path.\n";
        return ENOENT;
    }

    std::string filename{argv[1]};
    std::ifstream g_file(filename);
    if (!g_file.is_open()) {
        std::cout << "Can not open g-file.";
        return ENOENT;
    }

    GFileRawData g_file_data;
    g_file >> g_file_data;
    if (!g_file_data.is_complete()) { return EPERM; }
    g_file.close();

    auto t_after_read_file = high_resolution_clock::now();

    intp::InterpolationFunction<double, 2u> flux_function(
        3, g_file_data.flux,
        std::make_pair(g_file_data.r_center - .5 * g_file_data.dim.x(),
                       g_file_data.r_center + .5 * g_file_data.dim.x()),
        std::make_pair(g_file_data.z_mid - .5 * g_file_data.dim.y(),
                       g_file_data.z_mid + .5 * g_file_data.dim.y()));

    auto t_after_psi_mesh = high_resolution_clock::now();

    constexpr size_t RADIAL_GRID_COUNT = 129;
    constexpr size_t POLOIDAL_GRID_COUNT =
        255;  // need to be odd, required by the spline routine in GTC

    std::vector<Contour> contours;
    contours.reserve(RADIAL_GRID_COUNT);

    constexpr double PI2 = 2 * M_PI;

    double flux_full_range =
        g_file_data.flux_LCFS - g_file_data.flux_magnetic_axis;
    double flux_wall = .99 * flux_full_range;
    double flux_delta = flux_wall / (RADIAL_GRID_COUNT - 1);
    const double theta_delta = PI2 / POLOIDAL_GRID_COUNT;

    // Radial direction is divided into `RADIAL_GRID_COUNT-1` intervals.
    // `Contours` store `RADIAL_GRID_COUNT - 1` contours.
    for (size_t i = 0; i < RADIAL_GRID_COUNT; ++i) {
        const double psi = g_file_data.flux_magnetic_axis +
                           (i == 0 ? flux_delta
                            : i == RADIAL_GRID_COUNT - 1
                                ? flux_wall
                                : (static_cast<double>(i) + .5) * flux_delta);
        contours.emplace_back(psi, flux_function, g_file_data);
    }

    auto t_after_contour_construction = high_resolution_clock::now();

    // poloidal current
    intp::InterpolationFunction1D<> poloidal_current_intp{
        std::make_pair(0., flux_full_range),
        intp::util::get_range(g_file_data.f_pol), 3};

    // construct boozer coordinate, integrate B^2 * J along each contour
    auto b2j = [&](Vec<2, double> pt, double psi) {
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

    auto bp_field_square = [&](Vec<2, double> pt, double) {
        double dp_dr = flux_function.derivative(pt, {1, 0});
        double dp_dz = flux_function.derivative(pt, {0, 1});

        return (dp_dr * dp_dr + dp_dz * dp_dz) / (pt.x() * pt.x());
    };

    auto bt_field_square = [&](Vec<2, double> pt, double psi) {
        double f = poloidal_current_intp(psi);
        return f * f / (pt.x() * pt.x());
    };

    std::vector<double> poloidal_angles{g_file_data.geometric_poloidal_angles};
    poloidal_angles.push_back(poloidal_angles.front() + PI2);
    intp::InterpolationFunctionTemplate1D<> poloidal_template{
        intp::util::get_range(poloidal_angles), poloidal_angles.size(), 5,
        true};
    poloidal_angles.insert(poloidal_angles.begin(), 0);
    poloidal_angles.back() = PI2;
    intp::InterpolationFunctionTemplate1D<> poloidal_template_full{
        intp::util::get_range(poloidal_angles), poloidal_angles.size(), 5,
        false};

    // exclude magnetic axis; sample points are arranged such that knot points
    // coincide with output data points.

    intp::Mesh<double, 2> magnetic_boozer{RADIAL_GRID_COUNT,
                                          POLOIDAL_GRID_COUNT + 1};
    intp::Mesh<double, 2> r_boozer{RADIAL_GRID_COUNT, POLOIDAL_GRID_COUNT + 1};
    intp::Mesh<double, 2> z_boozer{RADIAL_GRID_COUNT, POLOIDAL_GRID_COUNT + 1};
    intp::Mesh<double, 2> jacobian_boozer{RADIAL_GRID_COUNT,
                                          POLOIDAL_GRID_COUNT + 1};

    std::vector<double> b2j_average;
    std::vector<double> normalized_toroidal_current;
    b2j_average.reserve(RADIAL_GRID_COUNT);
    normalized_toroidal_current.reserve(RADIAL_GRID_COUNT);

    // first point in radial direction is magnetic axis

    double B0 = b_field(g_file_data.magnetic_axis, 0.);
    double R0 = g_file_data.magnetic_axis.x();

#ifdef _DEBUG
    for (size_t ri = 0; ri < contours.size(); ++ri) {
        double psi = static_cast<double>(ri + 1) * flux_delta +
                     g_file_data.flux_magnetic_axis;
        // contour pt number + 1
        const size_t poloidal_size = contours[ri].size() + 1;
        std::vector<double> r_geo, z_geo, b2j_geo, bp_square_geo, bt_square_geo;
        std::vector<double> geometric_angles{
            g_file_data.geometric_poloidal_angles};
        geometric_angles.push_back(geometric_angles.front() + 2 * M_PI);
        std::vector<Contour::pt_type> contour_pts;
        contour_pts.reserve(poloidal_size);

        for (size_t i = 0; i < poloidal_size; ++i) {
            contour_pts.push_back(contours[ri][i % (poloidal_size - 1)]);
        }

        intp::InterpolationFunction1D<Contour::pt_type> contour_intp(
            intp::util::get_range(geometric_angles),
            intp::util::get_range(contour_pts), 5, true);
        auto psi_err_2 = [&](double t) {
            return std::pow(flux_function(contour_intp(t)) - psi, 2);
        };

        std::cout << psi << " -> "
                  << std::sqrt(util::integrate(psi_err_2, 0., 2 * M_PI) /
                               (2 * M_PI)) /
                         std::abs(psi)
                  << '\n';
    }
#endif

    std::vector<double> psi_sample_normalized;
    psi_sample_normalized.reserve(RADIAL_GRID_COUNT);
    // iterate through contour
    for (std::size_t ri = 0; ri < contours.size(); ++ri) {
        const double psi = contours[ri].flux() - g_file_data.flux_magnetic_axis;
        psi_sample_normalized.push_back(psi / (B0 * R0 * R0));
        // contour pt number + 1
        const size_t poloidal_size = contours[ri].size() + 1;
        std::vector<double> r_geo, z_geo, b2j_geo, bp_square_geo, bt_square_geo;
        r_geo.reserve(poloidal_size);
        z_geo.reserve(poloidal_size);
        b2j_geo.reserve(poloidal_size);
        bp_square_geo.reserve(poloidal_size);
        bt_square_geo.reserve(poloidal_size);
        for (size_t i = 0; i < poloidal_size; ++i) {
            const auto& pt = contours[ri][i % (poloidal_size - 1)];
            r_geo.push_back(pt.x());
            z_geo.push_back(pt.y());
            b2j_geo.push_back(b2j(pt, psi));
            bp_square_geo.push_back(bp_field_square(pt, psi) *
                                    j_field(pt, psi));
            bt_square_geo.push_back(bt_field_square(pt, psi) *
                                    j_field(pt, psi));
        }

        auto r_geo_intp =
            poloidal_template.interpolate(intp::util::get_range(r_geo));
        auto z_geo_intp =
            poloidal_template.interpolate(intp::util::get_range(z_geo));
        auto b2j_geo_intp =
            poloidal_template.interpolate(intp::util::get_range(b2j_geo));
        auto bp_square_geo_intp =
            poloidal_template.interpolate(intp::util::get_range(bp_square_geo));
        auto bt_square_geo_intp =
            poloidal_template.interpolate(intp::util::get_range(bt_square_geo));

        // integrate then interpolate to get \theta_boozer(\theta_geo)

        std::vector<double> boozer_geo;
        boozer_geo.reserve(poloidal_angles.size());
        boozer_geo.push_back(0);

        std::vector<double> bpa;
        bpa.reserve(poloidal_angles.size());
        bpa.push_back(0);
        std::vector<double> bta;
        bta.reserve(poloidal_angles.size());
        bta.push_back(0);
        // Poloidal grid begins from \theta = 0 and ends at \theta = 2\pi
        for (size_t i = 1; i < poloidal_angles.size(); ++i) {
            boozer_geo.push_back(boozer_geo.back() +
                                 util::integrate_coarse(b2j_geo_intp,
                                                        poloidal_angles[i - 1],
                                                        poloidal_angles[i]));
            bpa.push_back(bpa.back() +
                          util::integrate_coarse(bp_square_geo_intp,
                                                 poloidal_angles[i - 1],
                                                 poloidal_angles[i]));
            bta.push_back(bta.back() +
                          util::integrate_coarse(bt_square_geo_intp,
                                                 poloidal_angles[i - 1],
                                                 poloidal_angles[i]));
        }
        b2j_average.push_back(boozer_geo.back() / PI2);
        normalized_toroidal_current.push_back(bpa.back() / (PI2 * R0 * B0));
        // normalization
        for (auto& v : boozer_geo) { v /= b2j_average.back(); }
        auto boozer_geo_intp = poloidal_template_full.interpolate(
            intp::util::get_range(boozer_geo));

        // calculate necessary values on a even-spaced boozer grid
        for (size_t i = 0; i <= POLOIDAL_GRID_COUNT; ++i) {
            double theta_boozer =
                (static_cast<double>(i % POLOIDAL_GRID_COUNT) + .5) *
                theta_delta;
            double theta_geo = util::find_root(
                [&](double t) { return boozer_geo_intp(t) - theta_boozer; }, 0.,
                PI2);
            double r_grid = r_geo_intp(theta_geo);
            double z_grid = z_geo_intp(theta_geo);

            // be careful of normalization

            magnetic_boozer(ri, i) = b_field({r_grid, z_grid}, psi) / B0;
            r_boozer(ri, i) = r_grid / R0;
            // z value is shifted such that magnetic axis has z = 0
            z_boozer(ri, i) = (z_grid - g_file_data.magnetic_axis.y()) / R0;
            jacobian_boozer(ri, i) = j_field({r_grid, z_grid}, psi) / R0 * B0;
        }
    }

    auto t_after_boozer_mesh = high_resolution_clock::now();

    // from now on, magnetic flux should be normalized

    flux_full_range /= B0 * R0 * R0;
    flux_wall /= B0 * R0 * R0;
    flux_delta /= B0 * R0 * R0;

    // interpolation of 2D functions on boozer grid

#define create_2d_spline(field)                                 \
    intp::InterpolationFunction<double, 2> field##_boozer_intp( \
        2, {false, true}, field##_boozer,                       \
        intp::util::get_range(psi_sample_normalized),           \
        std::make_pair(.5 * theta_delta, PI2 + .5 * theta_delta))

    create_2d_spline(magnetic);
    create_2d_spline(r);
    create_2d_spline(z);
    create_2d_spline(jacobian);

    auto t_after_interpolation_on_boozer_mesh = high_resolution_clock::now();

    // interpolation of 1D functions

    intp::InterpolationFunction1D<> safety_factor_intp(
        std::make_pair(0., flux_full_range),
        intp::util::get_range(g_file_data.safety_factor), 3);

    // resample safety factor
    std::vector<double> safety_factor_resample(RADIAL_GRID_COUNT);
    for (std::size_t i = 0; i < RADIAL_GRID_COUNT; ++i) {
        // safety_factor_resample[i] =
        //     safety_factor_intp(static_cast<double>(i) * flux_delta);
        safety_factor_resample[i] =
            safety_factor_intp(psi_sample_normalized[i]);
    }

    intp::InterpolationFunction1D<> safety_factor_resample_intp(
        intp::util::get_range(psi_sample_normalized),
        intp::util::get_range(safety_factor_resample), 2);

    // resample poloidal current
    std::vector<double> normalized_poloidal_current_resample(RADIAL_GRID_COUNT);
    for (std::size_t i = 0; i < RADIAL_GRID_COUNT; ++i) {
        normalized_poloidal_current_resample[i] =
            poloidal_current_intp(psi_sample_normalized[i] * B0 * R0 * R0) /
            (B0 * R0);
    }

    intp::InterpolationFunction1D<> normalized_poloidal_current_resample_intp(
        intp::util::get_range(psi_sample_normalized),
        intp::util::get_range(normalized_poloidal_current_resample), 2);

    intp::InterpolationFunction1D<> normalized_toroidal_current_intp(
        intp::util::get_range(psi_sample_normalized),
        intp::util::get_range(normalized_toroidal_current), 2);

    constexpr double magnetic_constant = 4.e-7 * M_PI;
    intp::InterpolationFunction1D<> pressure_intp(
        std::make_pair(0., flux_full_range),
        intp::util::get_range(g_file_data.pressure), 3);

    // resample pressure profile
    std::vector<double> normalized_pressure_resample(RADIAL_GRID_COUNT);
    for (std::size_t i = 0; i < RADIAL_GRID_COUNT; ++i) {
        normalized_pressure_resample[i] =
            pressure_intp(psi_sample_normalized[i]) /
            (B0 * B0 / magnetic_constant);
    }

    intp::InterpolationFunction1D<> normalized_pressure_resample_intp(
        intp::util::get_range(psi_sample_normalized),
        intp::util::get_range(normalized_pressure_resample), 2);

    // integrate safety factor wrt poloidal flux to obtain toroidal flux

    std::vector<double> normalized_toroidal_flux;
    normalized_toroidal_flux.reserve(RADIAL_GRID_COUNT);
    for (size_t i = 0; i < RADIAL_GRID_COUNT; ++i) {
        normalized_toroidal_flux.push_back(
            (i == 0 ? 0. : normalized_toroidal_flux.back()) +
            util::integrate_coarse(safety_factor_intp,
                                   i == 0 ? 0. : psi_sample_normalized[i - 1],
                                   psi_sample_normalized[i]));
    }
    intp::InterpolationFunction1D<> normalized_toroidal_flux_intp(
        intp::util::get_range(psi_sample_normalized),
        intp::util::get_range(normalized_toroidal_flux), 2);

    std::vector<double> normalized_r_minor(std::move(normalized_toroidal_flux));
    for (auto& v : normalized_r_minor) {
        // define r as if \psi_t = B_0 * r^2 / 2
        v = std::sqrt(2. * v);
    }
    intp::InterpolationFunction1D<> normalized_r_minor_intp(
        intp::util::get_range(psi_sample_normalized),
        intp::util::get_range(normalized_r_minor), 2);

    auto t_after_1d_interpolations = high_resolution_clock::now();

    // consistency check

    // define output functions for 1D and 2D spline coefficients

    auto write_1d_coef = [&](std::ostream& os, const auto& f_1d,
                             std::size_t idx, double value_on_axis = 0.,
                             bool singular = false) {
        // const cp =
        double c0, c1, c2;
        if (idx == 0) {
            const double v1 = f_1d(flux_delta) - value_on_axis;
            const double v2 = f_1d.derivative(std::make_pair(flux_delta, 1)) *
                              (singular ? 2. * std::sqrt(flux_delta) : 1.);
            const double d = singular ? std::sqrt(flux_delta) : flux_delta;
            c0 = value_on_axis;
            c1 = 2. * v1 / d - v2;
            c2 = (-v1 + v2 * d) / (d * d);
        } else {
            const double psi = static_cast<double>(idx) * flux_delta;
            c0 = f_1d(psi);
            c1 = f_1d.derivative(std::make_pair(psi, 1));
            c2 = .5 * f_1d.derivative(std::make_pair(psi + .5 * flux_delta, 2));
        }
        os << std::setw(18) << c0 << std::setw(18) << c1 << std::setw(18) << c2
           << '\n';
    };

    auto write_2d_coef = [&](std::ostream& os, auto& f_2d, std::size_t idx,
                             double value_on_axis = 0.) {
        for (size_t i = 0; i < 9; ++i) {
            size_t psi_order = i % 3;
            size_t theta_order = i / 3;
            for (size_t j = 0; j <= POLOIDAL_GRID_COUNT; ++j) {
                double coef =
                    (psi_order == 2 ? .5 : 1.) * (theta_order == 2 ? .5 : 1.);
                const double theta =
                    (static_cast<double>(j) + (theta_order == 2 ? .5 : 0.)) *
                    theta_delta;
                if (idx == 0) {
                    const double v1 =
                        theta_order == 0
                            ? f_2d(flux_delta, theta) - value_on_axis
                            : f_2d.derivative({flux_delta, theta},
                                              {0, theta_order});
                    const double v2 =
                        f_2d.derivative({flux_delta, theta}, {1, theta_order});

                    coef = psi_order == 0
                               ? theta_order == 0 ? value_on_axis : 0.
                               : 2. * coef *
                                     (psi_order == 1
                                          ? (v1 - v2 * flux_delta) /
                                                std::sqrt(flux_delta)
                                          : -(v1 - 2. * v2 * flux_delta) /
                                                flux_delta);
                } else {
                    const double psi = (static_cast<double>(idx) +
                                        (psi_order == 2 ? .5 : 0.)) *
                                       flux_delta;
                    coef = (i == 0 ? f_2d(psi, theta)
                                   : coef * f_2d.derivative(
                                                {psi, theta},
                                                {psi_order, theta_order}));
                }

                os << std::setw(18) << coef;
                if (j % 4 == 3) { os << '\n'; }
            }
            if ((POLOIDAL_GRID_COUNT + 1) % 4 != 0) { os << '\n'; }
        }
    };

    {
        std::ofstream sp_data("./spdata.dat", std::ios::out);
        sp_data << "Generated by MCT, from " << filename << '\n';
        sp_data << std::setw(4) << RADIAL_GRID_COUNT << std::setw(4)
                << POLOIDAL_GRID_COUNT << std::setw(4) << 4 << std::setw(4) << 8
                << '\n';
        sp_data << std::scientific << std::uppercase << std::setprecision(10);
        sp_data << std::setw(18) << flux_wall << std::setw(18) << flux_wall
                << '\n';

        const double g0n = g_file_data.f_pol[0] / (B0 * R0);
        const double q0 = g_file_data.safety_factor[0];
        const double p0n =
            g_file_data.pressure[0] / (B0 * B0 / magnetic_constant);
        for (size_t i = 0; i < RADIAL_GRID_COUNT; ++i) {
            write_2d_coef(sp_data, magnetic_boozer_intp, i, 1.);
            write_2d_coef(sp_data, r_boozer_intp, i, 1.);
            write_2d_coef(sp_data, z_boozer_intp, i);
            write_2d_coef(sp_data, jacobian_boozer_intp, i, g0n * q0);
            write_1d_coef(sp_data, safety_factor_resample_intp, i, q0);
            write_1d_coef(sp_data, normalized_poloidal_current_resample_intp, i,
                          g0n);
            write_1d_coef(sp_data, normalized_toroidal_current_intp, i);
            write_1d_coef(sp_data, normalized_pressure_resample_intp, i, p0n);
            write_1d_coef(sp_data, normalized_r_minor_intp, i, 0., true);
            write_1d_coef(sp_data, normalized_toroidal_flux_intp, i);
        }

        // TODO: ripple related
        sp_data << std::setw(4) << 0 << std::setw(4) << 0 << '\n';
        sp_data << std::setw(18) << R0 << std::setw(18) << 0. << std::setw(18)
                << 0. << '\n';
        sp_data << std::setw(18) << 0. << std::setw(18) << 0. << '\n';
    }

    auto t_after_output = high_resolution_clock::now();

    std::cout << "Creating spdata from g_file finished.\n";
    std::cout << "R0 = " << R0 << "m, B0 = " << B0 << "T\n";

    std::cout << "\nPhase\t\t\tTime consumption (ms)\n";
    std::cout << "read file\t\t"
              << duration<double, milliseconds::period>(t_after_read_file -
                                                        t_start)
                     .count()
              << '\n';
    std::cout << "psi interp\t\t"
              << duration<double, milliseconds::period>(t_after_psi_mesh -
                                                        t_after_read_file)
                     .count()
              << '\n';
    std::cout << "contours  \t\t"
              << duration<double, milliseconds::period>(
                     t_after_contour_construction - t_after_psi_mesh)
                     .count()
              << '\n';
    std::cout << "boozer mesh  \t\t"
              << duration<double, milliseconds::period>(
                     t_after_boozer_mesh - t_after_contour_construction)
                     .count()
              << '\n';
    std::cout << "2d interp\t\t"
              << duration<double, milliseconds::period>(
                     t_after_interpolation_on_boozer_mesh - t_after_boozer_mesh)
                     .count()
              << '\n';
    std::cout << "1d interp\t\t"
              << duration<double, milliseconds::period>(
                     t_after_1d_interpolations -
                     t_after_interpolation_on_boozer_mesh)
                     .count()
              << '\n';
    std::cout << "output   \t\t"
              << duration<double, milliseconds::period>(
                     t_after_output - t_after_1d_interpolations)
                     .count()
              << '\n';

    std::cout
        << "\nTotal time consumption: "
        << duration<double, seconds::period>(t_after_output - t_start).count()
        << "s\n";

    return 0;
}
