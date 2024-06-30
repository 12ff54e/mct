#include <cerrno>
#include <filesystem>
#include <fstream>
#include <iomanip>

#include "include/Contour.hpp"
#include "include/GFileRawData.hpp"
#include "include/Spdata.hpp"

#define ZQ_TIMER_IMPLEMENTATION
#include "include/Timer.h"

#include "include/clap.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Input {
    bool use_SI;
    unsigned radial_grid_count;
    unsigned poloidal_grid_count;
    std::string input_path;
    std::string output_path;
};

int main(int argc, char** argv) {
    auto& timer = Timer::get_timer();
    timer.start("Read g-file");

    CLAP_BEGIN(Input)
    CLAP_REGISTER(input_path, "--input-path")
    CLAP_REGISTER(output_path, "--output-path")
    CLAP_REGISTER(radial_grid_count, "--radial")
    CLAP_REGISTER(poloidal_grid_count, "--poloidal")
    CLAP_REGISTER(use_SI, "--use-si")
    CLAP_END(Input)

    Input config;
    config.use_SI = false;
    config.radial_grid_count = 129;
    config.poloidal_grid_count = 255;
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

    Spdata spdata(g_file_data, config.radial_grid_count,
                  config.poloidal_grid_count, config.use_SI);

    timer.pause_last_and_start_next("Write to spdata ");

    std::ofstream output(config.output_path, std::ios::out);
    output << "Generated by MCT, from " << filename << '\n';
    output << spdata;

    timer.pause();
    timer.print();

    return 0;
}
