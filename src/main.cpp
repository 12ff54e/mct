#include <cerrno>
#include <filesystem>
#include <fstream>
#include <iomanip>

#include "lib/Contour.h"
#include "lib/GFileRawData.h"
#include "lib/Spdata.h"

#define ZQ_TIMER_IMPLEMENTATION
#include "lib/Timer.h"

#include "lib/clap.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Input {
    bool use_SI;
    unsigned radial_grid_count;
    unsigned poloidal_grid_count;
    double psi_ratio;
    std::string input_path;
    std::string output_path;
};

int main(int argc, char** argv) {
    auto& timer = Timer::get_timer();
    timer.start("Read g-file");

    CLAP_BEGIN(Input)
    CLAP_ADD_USAGE("[OPTION]... INPUT_FILE")
    CLAP_ADD_DESCRIPTION("Transfrom g-file to spdata using in GTC.")
    CLAP_REGISTER_ARG(input_path)
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(
        output_path, "--output-path", "-o",
        "specify the path of output file, default to '$PWD/spdata.dat'")
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(radial_grid_count, "--radial", "-r",
                                          "radial grid number")
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(poloidal_grid_count, "--poloidal",
                                          "-p", "poloidal grid number")
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(
        use_SI, "--use-si",
        "output with SI unit instead of normalized to B_0 and R_0")
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(
        psi_ratio, "--psi-ratio",
        "ratio of psiw in output to that given in gfile, default to be 0.99")
    CLAP_END(Input)

    Input config;
    config.use_SI = false;
    config.radial_grid_count = 129;
    config.poloidal_grid_count = 255;
    config.psi_ratio = 0.99;
    try {
        CLAP<Input>::parse_input(config, argc, argv);
    } catch (std::exception& e) {
        std::cerr << e.what();
        return EINVAL;
    }

    if (config.input_path.empty()) {
        std::cerr << "Please provide input file.";
        return ENOENT;
    }
    if (config.output_path.empty()) {
        config.output_path =
            (std::filesystem::current_path() / "spdata.dat").string();
    }

    const auto& filename = config.input_path;
    std::ifstream g_file(filename);
    if (!g_file.is_open()) {
        std::cerr << "Can not open g-file.\n";
        return ENOENT;
    }

    GFileRawData g_file_data;
    g_file >> g_file_data;
    if (!g_file_data.is_complete()) { return EPERM; }
    g_file.close();

    timer.pause_last_and_start_next("Generate BSpline");

    constexpr std::size_t radial_sample_point = 256;
    Spdata spdata(g_file_data, config.radial_grid_count,
                  config.poloidal_grid_count, config.use_SI,
                  radial_sample_point, config.psi_ratio);

    timer.pause_last_and_start_next("Write to spdata ");

    std::ofstream output(config.output_path, std::ios::out);
    output << "Generated by MCT, from " << filename << '\n';
    output << spdata;

    double R0 = g_file_data.magnetic_axis.x();
    std::cout << "R0 = " << R0
              << "m\nB0 = " << std::abs(g_file_data.f_pol[0] / R0) << "T\n";

    timer.pause();
    timer.print();

    return 0;
}
