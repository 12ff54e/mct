#include <cmath>
#include <iostream>
#include <limits>

#include "../src/include/util.hpp"
#include "Assertion.hpp"

int main(int argc, char const* argv[]) {
    Assertion assertion;
    assertion(
        util::find_root([](double x) { return x * x - 2; }, 1., 2.) / M_SQRT2 -
                1. <=
            std::numeric_limits<double>::epsilon(),
        "find root function not correct!");
    return assertion.status();
}
