#include <chrono>
#include <fstream>

#include "../src/include/Contour.hpp"
#include "../src/include/GFileRawData.hpp"
#include "Assertion.hpp"

int main(int argc, char const* argv[]) {
    using namespace std::chrono;
    Assertion assertion;

    auto t_start = high_resolution_clock::now();

    std::ifstream g_file;
    g_file.open("../data/gfile-cfetr5.7-baseline-129", std::ios::in);
    if (!g_file.is_open()) {
        assertion(false, "Can not open g-file.");
        return assertion.status();
    }

    GFileRawData g_file_data;
    g_file >> g_file_data;
    assertion(g_file_data.is_complete(), "Parse g-file failed.");
    if (assertion.last_status()) { return assertion.status(); }

    auto t_after_read_file = high_resolution_clock::now();

    intp::InterpolationFunction<double, 2u> flux_function(
        3, g_file_data.flux,
        std::make_pair(g_file_data.r_center - .5 * g_file_data.dim.x(),
                       g_file_data.r_center + .5 * g_file_data.dim.x()),
        std::make_pair(g_file_data.z_mid - .5 * g_file_data.dim.y(),
                       g_file_data.z_mid + .5 * g_file_data.dim.y()));

    auto t_after_psi_mesh = high_resolution_clock::now();

    double psi = .5 * (g_file_data.flux_magnetic_axis + g_file_data.flux_LCFS);
    Contour middle_contour(psi, flux_function, g_file_data);

    auto t_after_middle_contour_construction = high_resolution_clock::now();

    double err{};
    for (size_t i = 0; i < middle_contour.size(); ++i) {
        double local_err =
            std::abs(flux_function(middle_contour[i]) / psi - 1.);
        err += local_err;
#ifdef _DEBUG
        std::cout << "[DEBUG] " << local_err << '\n';
#endif
    }
    err /= middle_contour.size();

    assertion(err < 1e-15, "Middle contour deviation is too much.\n");

    std::cout << "Contour test "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative error ~ " << err << '\n';

    // If the integrand depends on geometric theta only, the result should be
    // the same for any convex contour around magnetic axis.
    auto one_minus_cos = middle_contour.indefinite_integrate_along(
        [&](Vec<2, double> pt) {
            return std::sin(util::arctan(pt - g_file_data.magnetic_axis));
        },
        false);

    err = 0;
    constexpr size_t sample_num = 100;
    for (size_t i = 0; i < sample_num; ++i) {
        double theta = 2 * M_PI * i / sample_num;
        double local_err =
            std::abs(one_minus_cos(theta) - 1. + std::cos(theta));
        err += local_err;
#ifdef _DEBUG
        std::cout << "[DEBUG] " << local_err << '\n';
#endif
    }
    err /= sample_num;

    assertion(err < 1e-6,
              "Integration along contour did not work as expected.\n");

    std::cout << "Contour integration test "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative error ~ " << err << '\n';

    std::cout << "Phase\t\t\tTime consumption (ms)\n";
    std::cout
        << "read file\t\t"
        << duration_cast<microseconds>(t_after_read_file - t_start).count() /
               1000.
        << '\n';
    std::cout << "psi interp\t\t"
              << duration_cast<microseconds>(t_after_psi_mesh -
                                             t_after_read_file)
                         .count() /
                     1000.
              << '\n';
    std::cout << "contour  \t\t"
              << duration_cast<microseconds>(
                     t_after_middle_contour_construction - t_after_psi_mesh)
                         .count() /
                     1000.
              << '\n';

    return assertion.status();
}
