#include <fstream>

#include "../src/include/GFileRawData.hpp"
#include "Assertion.hpp"

int main(int argc, char const* argv[]) {
    Assertion assertion;

    std::ifstream g_file;
    g_file.open("../data/gfile-cfetr5.7-baseline-129", std::ios::in);
    if (g_file.is_open()) {
        GFileRawData g_file_data;
        g_file >> g_file_data;
        assertion(g_file_data.is_complete(), "Parse gfile failed.");
    } else {
        assertion(false, "File open failed.");
    }

    return assertion.status();
}
