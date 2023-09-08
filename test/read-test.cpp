#include <chrono>
#include <fstream>
#include <iomanip>

#include "../src/include/Contour.hpp"
#include "../src/include/GFileRawData.hpp"
#include "Assertion.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    using namespace std::chrono;
    Assertion assertion;

    auto t_start = high_resolution_clock::now();

    std::string filename = "gfile-cfetr5.7-baseline-129";
    std::ifstream g_file(std::string("../data/") + filename);
    if (!g_file.is_open()) {
        assertion(false, "Can not open g-file.");
        return assertion.status();
    }

    GFileRawData g_file_data;
    g_file >> g_file_data;
    assertion(g_file_data.is_complete(), "Parse g_file failed.");
    if (assertion.last_status()) { return assertion.status(); }
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
    constexpr size_t POLOIDAL_GRID_COUNT = 255;

    std::vector<Contour> contours;
    contours.reserve(RADIAL_GRID_COUNT - 1);

    constexpr double PI2 = 2 * M_PI;

    double flux_diff = g_file_data.flux_LCFS - g_file_data.flux_magnetic_axis;
    double flux_wall = .99 * flux_diff;
    double flux_delta = flux_wall / (RADIAL_GRID_COUNT - 1);
    const double theta_delta = PI2 / POLOIDAL_GRID_COUNT;

    // Radial direction is divided into `RADIAL_GRID_COUNT-1` intervals.
    // `Contours` store `RADIAL_GRID_COUNT - 1` contours.
    for (size_t i = 1; i < RADIAL_GRID_COUNT; ++i) {
        double psi = g_file_data.flux_magnetic_axis + i * flux_delta;
        contours.emplace_back(psi, flux_function, g_file_data);
    }

    auto t_after_contour_construction = high_resolution_clock::now();

    // poloidal current
    intp::InterpolationFunction1D<> poloidal_current_intp{
        std::make_pair(0., flux_diff), intp::util::get_range(g_file_data.f_pol),
        3};

    // construct boozer coordinate, integrate B^2 * J along each contour
    auto b2j = [&](Vec<2, double> pt, double psi) {
        double dp_dr = flux_function.derivative(pt, {1, 0});
        double dp_dz = flux_function.derivative(pt, {0, 1});
        double r = pt.x();
        pt -= g_file_data.magnetic_axis;
        double r2 = pt.__L2_norm_square();
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
        double r2 = pt.__L2_norm_square();

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

    intp::Mesh<double, 2> magnetic_boozer{RADIAL_GRID_COUNT,
                                          POLOIDAL_GRID_COUNT + 1};
    intp::Mesh<double, 2> r_boozer{RADIAL_GRID_COUNT, POLOIDAL_GRID_COUNT + 1};
    intp::Mesh<double, 2> z_boozer{RADIAL_GRID_COUNT, POLOIDAL_GRID_COUNT + 1};
    intp::Mesh<double, 2> jacobian_boozer{RADIAL_GRID_COUNT,
                                          POLOIDAL_GRID_COUNT + 1};

    std::vector<double> b2j_average;
    std::vector<double> tc;
    b2j_average.reserve(RADIAL_GRID_COUNT);
    tc.reserve(RADIAL_GRID_COUNT);

    // first point in radial direction is magnetic axis

    double B0 = b_field(g_file_data.magnetic_axis, 0.);
    double R0 = g_file_data.magnetic_axis.x();
    {
        auto set_ma_val = [&](double val, auto& field) {
            for (size_t i = 0; i < POLOIDAL_GRID_COUNT; ++i) {
                field(0, i) = val;
            }
        };
        set_ma_val(1., magnetic_boozer);
        set_ma_val(1., r_boozer);
        set_ma_val(0., z_boozer);
        set_ma_val(0., jacobian_boozer);

        b2j_average.push_back(g_file_data.safety_factor[0] *
                              g_file_data.f_pol[0]);
        tc.push_back(0.);
    }
    // iterate through contour
    for (size_t ri = 0; ri < contours.size(); ++ri) {
        double psi = (ri + 1) * flux_delta;
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
        tc.push_back(bpa.back() / (PI2 * R0 * B0));
        // normalization
        for (auto& v : boozer_geo) { v /= b2j_average.back(); }
        auto boozer_geo_intp = poloidal_template_full.interpolate(
            intp::util::get_range(boozer_geo));

        // calculate necessary values on a even-spaced boozer grid
        for (size_t i = 0; i <= POLOIDAL_GRID_COUNT; ++i) {
            double theta_boozer = i * theta_delta;
            double theta_geo =
                i == 0 || i == POLOIDAL_GRID_COUNT
                    ? theta_boozer
                    : util::find_root(
                          [&](double t) {
                              return boozer_geo_intp(t) - theta_boozer;
                          },
                          0., PI2);
            double r_grid = r_geo_intp(theta_geo);
            double z_grid = z_geo_intp(theta_geo);

            // be careful of normalization

            magnetic_boozer(ri + 1, i) = b_field({r_grid, z_grid}, psi) / B0;
            r_boozer(ri + 1, i) = r_grid / R0;
            // z value is shifted such that magnetic axis has z = 0
            z_boozer(ri + 1, i) = (z_grid - g_file_data.magnetic_axis.y()) / R0;
            jacobian_boozer(ri + 1, i) =
                j_field({r_grid, z_grid}, psi) / R0 * B0;
        }
    }

    auto t_after_boozer_mesh = high_resolution_clock::now();

    // from now on, magnetic flux should be normalized

    flux_diff /= B0 * R0 * R0;
    flux_wall /= B0 * R0 * R0;
    flux_delta /= B0 * R0 * R0;

    // interpolation of 2D functions on boozer grid

#define create_2d_spline(field)                                          \
    intp::InterpolationFunction<double, 2> field##_boozer_intp(          \
        2, {false, true}, field##_boozer, std::make_pair(0., flux_wall), \
        std::make_pair(0., PI2))

    create_2d_spline(magnetic);
    create_2d_spline(r);
    create_2d_spline(z);
    create_2d_spline(jacobian);

    auto t_after_interpolation_on_boozer_mesh = high_resolution_clock::now();

    // interpolation of 1D functions

    intp::InterpolationFunction1D<> safety_factor_intp(
        std::make_pair(0., flux_wall),
        intp::util::get_range(g_file_data.safety_factor), 2);

    auto normalized_poloidal_current = g_file_data.f_pol;
    for (auto& v : normalized_poloidal_current) { v /= (B0 * R0); }
    intp::InterpolationFunction1D<> normalized_poloidal_current_intp(
        std::make_pair(0., flux_diff),
        intp::util::get_range(normalized_poloidal_current), 2);

    std::vector<double> normalized_toroidal_current(b2j_average);
    for (size_t i = 0; i < normalized_toroidal_current.size(); ++i) {
        normalized_toroidal_current[i] =
            normalized_toroidal_current[i] / (B0 * R0) -
            safety_factor_intp(i * flux_delta) *
                normalized_poloidal_current_intp(i * flux_delta);
    }
    intp::InterpolationFunction1D<> normalized_toroidal_current_intp(
        std::make_pair(0., flux_wall), intp::util::get_range(tc), 2);

    constexpr double magnetic_constant = 4.e-7 * M_PI;
    auto normalized_pressure = g_file_data.pressure;
    for (auto& v : normalized_pressure) { v /= (B0 * B0 / magnetic_constant); }
    intp::InterpolationFunction1D<> normalized_pressure_intp(
        std::make_pair(0., flux_diff),
        intp::util::get_range(normalized_pressure), 2);

    // integrate safety factor wrt poloidal flux to obtain toroidal flux

    std::vector<double> normalized_toroidal_flux;
    normalized_toroidal_flux.reserve(RADIAL_GRID_COUNT);
    normalized_toroidal_flux.push_back(0);
    for (size_t i = 1; i < RADIAL_GRID_COUNT; ++i) {
        normalized_toroidal_flux.push_back(
            normalized_toroidal_flux.back() +
            util::integrate_coarse(safety_factor_intp, flux_delta * (i - 1),
                                   flux_delta * i));
    }
    intp::InterpolationFunction1D<> normalized_toroidal_flux_intp(
        std::make_pair(0., flux_wall),
        intp::util::get_range(normalized_toroidal_flux), 2);

    std::vector<double> normalized_r_minor(std::move(normalized_toroidal_flux));
    for (auto& v : normalized_r_minor) {
        // define r as if \psi_t = B_0 * r^2 / 2
        v = std::sqrt(2. * v);
    }
    intp::InterpolationFunction1D<> normalized_r_minor_intp(
        std::make_pair(0., flux_wall),
        intp::util::get_range(normalized_r_minor), 2);

    auto t_after_1d_interpolations = high_resolution_clock::now();

    // define output functions for 1D and 2D spline coefficients

    auto write_1d_coef = [](std::ostream& os, auto& f_1d, double psi) {
        os << std::setw(18) << f_1d(psi) << std::setw(18)
           << f_1d.derivative(std::make_pair(psi, 1)) << std::setw(18)
           << .5 * f_1d.derivative(std::make_pair(psi, 2)) << '\n';
    };

    auto write_2d_coef = [&](std::ostream& os, auto& f_2d, double psi) {
        for (size_t i = 0; i < 9; ++i) {
            size_t psi_order = i % 3;
            size_t theta_order = i / 3;
            double coef =
                (psi_order == 2 ? .5 : 1.) * (theta_order == 2 ? .5 : 1.);
            for (size_t j = 0; j <= POLOIDAL_GRID_COUNT; ++j) {
                os << std::setw(18)
                   << (i == 0
                           ? f_2d(psi, theta_delta * j)
                           : coef * f_2d.derivative({psi, j * theta_delta},
                                                    {psi_order, theta_order}));
                if (j % 4 == 3) { os << '\n'; }
            }
            if ((POLOIDAL_GRID_COUNT + 1) % 4 != 0) { os << '\n'; }
        }
    };
    auto write_2d_coef_first_seg = [&](std::ostream& os, auto& f_2d) {
        for (size_t i = 0; i < 9; ++i) {
            size_t psi_order = i % 3;
            size_t theta_order = i / 3;
            for (size_t j = 0; j <= POLOIDAL_GRID_COUNT; ++j) {
                double v1 = theta_order == 0
                                ? f_2d(flux_delta, theta_delta * j) -
                                      f_2d(0., theta_delta * j)
                                : f_2d.derivative({flux_delta, theta_delta * j},
                                                  {0, theta_order}) -
                                      f_2d.derivative({0., theta_delta * j},
                                                      {0, theta_order});
                double v2 = f_2d.derivative({flux_delta, theta_delta * j},
                                            {1, theta_order});
                std::array<double, 3> psi_coef{
                    theta_order == 0 ? f_2d(0., theta_delta * j) : 0.,
                    2. * (v1 - v2 * flux_delta) / std::sqrt(flux_delta),
                    -(v1 - 2. * v2 * flux_delta) / flux_delta};
                os << std::setw(18)
                   << (theta_order == 2 ? .5 : 1.) * psi_coef[psi_order];
                if (j % 4 == 3) { os << '\n'; }
            }
            if ((POLOIDAL_GRID_COUNT + 1) % 4 != 0) { os << '\n'; }
        }
    };

    if constexpr (true) {
        std::ofstream sp_data("../data/spdata.dat", std::ios::out);
        sp_data << "Generated by MCT, from " << filename << '\n';
        sp_data << std::setw(4) << RADIAL_GRID_COUNT << std::setw(4)
                << POLOIDAL_GRID_COUNT << std::setw(4) << 4 << std::setw(4) << 8
                << '\n';
        sp_data << std::scientific << std::uppercase << std::setprecision(10);
        sp_data << std::setw(18) << flux_wall << std::setw(18) << flux_wall
                << '\n';
        for (size_t i = 0; i < RADIAL_GRID_COUNT; ++i) {
            // first interval is different
            if (i == 0) {
                write_2d_coef_first_seg(sp_data, magnetic_boozer_intp);
                write_2d_coef_first_seg(sp_data, r_boozer_intp);
                write_2d_coef_first_seg(sp_data, z_boozer_intp);
                write_2d_coef_first_seg(sp_data, jacobian_boozer_intp);
            } else {
                write_2d_coef(sp_data, magnetic_boozer_intp, i * flux_delta);
                write_2d_coef(sp_data, r_boozer_intp, i * flux_delta);
                write_2d_coef(sp_data, z_boozer_intp, i * flux_delta);
                write_2d_coef(sp_data, jacobian_boozer_intp, i * flux_delta);
            }
            write_1d_coef(sp_data, safety_factor_intp, i * flux_delta);
            write_1d_coef(sp_data, normalized_poloidal_current_intp,
                          i * flux_delta);
            write_1d_coef(sp_data, normalized_toroidal_current_intp,
                          i * flux_delta);
            write_1d_coef(sp_data, normalized_pressure_intp, i * flux_delta);
            // r_minor has different definition of coefficient at first interval
            if (i == 0) {
                double v1 = normalized_r_minor_intp(flux_delta);
                double v2 = normalized_r_minor_intp.derivative({flux_delta}, 1);
                double c1 = 2. * (v1 - v2 * flux_delta) / std::sqrt(flux_delta);
                double c2 = -(v1 - 2. * v2 * flux_delta) / flux_delta;
                sp_data << std::setw(18) << 0. << std::setw(18) << c1
                        << std::setw(18) << c2 << '\n';
            } else {
                write_1d_coef(sp_data, normalized_r_minor_intp, i * flux_delta);
            }
            write_1d_coef(sp_data, normalized_toroidal_flux_intp,
                          i * flux_delta);
        }
        // TODO: ripple related
        sp_data << std::setw(4) << 0 << std::setw(4) << 0 << '\n';  // d0, brip
        sp_data << std::setw(18) << R0 << std::setw(18) << 0. << std::setw(18)
                << 0. << '\n';  // R_major, d0, brip
        sp_data << std::setw(18) << 0. << std::setw(18) << 0.
                << '\n';  // wrip, xrip
    }

    auto t_after_output = high_resolution_clock::now();

    std::cout << "Creating spdata from g_file finished.\n";
    std::cout << "R0 = " << R0 << ", B0 = " << B0 << '\n';

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

    std::cout << B0 * B0 * R0 /
                     std::sqrt(flux_function.derivative(
                                   g_file_data.magnetic_axis, {2, 0}) *
                               flux_function.derivative(
                                   g_file_data.magnetic_axis, {0, 2}))
              << ", " << g_file_data.safety_factor[0] * g_file_data.f_pol[0]
              << '\n';

    return assertion.status();
}
