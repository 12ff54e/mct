#include <algorithm>
#include <iostream>

#include "include/GFileRawData.hpp"
#include "include/util.hpp"

bool GFileRawData::is_complete() const noexcept {
    return __complete;
}

std::ifstream& operator>>(std::ifstream& is, GFileRawData& g) {
    std::string s_dum;
    double d_dum;
    // 1st line
    is >> g.code_name >> g.date >> s_dum;
    is >> g.shot_id >> g.time >> s_dum;
    if (is.tellg() < 51) { is >> s_dum; }
    is >> g.nw >> g.nh;
    if (is.fail()) {
        std::cout << "Data corruption at 1st line.\n";
        return is;
    }
    // 2nd line
    is >> g.dim.x() >> g.dim.y() >> g.r_center >> g.r_left >> g.z_mid;
    if (is.fail()) {
        std::cout << "Data corruption at 2nd line.\n";
        return is;
    }
    // 3rd line
    is >> g.magnetic_axis.x() >> g.magnetic_axis.y() >> g.flux_magnetic_axis >>
        g.flux_LCFS >> g.b_center;
    if (is.fail()) {
        std::cout << "Data corruption at 3rd line.\n";
        return is;
    }
    // 4th line
    is >> g.current >> d_dum >> d_dum >> d_dum >> d_dum;
    if (is.fail()) {
        std::cout << "Data corruption at 4th line.\n";
        return is;
    }
    // 5th line
    is >> d_dum >> d_dum >> g.flux_sep >> g.sep.x() >> g.sep.y();
    if (is.fail()) {
        std::cout << "Data corruption at 5th line.\n";
        return is;
    }

    auto read_vec = [&is](unsigned int count, std::vector<double>& arr) {
        double d;
        for (unsigned i = 0; i < count; ++i) {
            is >> d;
            arr.emplace_back(d);
        }
    };

    // poloidal current
    read_vec(g.nw, g.f_pol);
    if (is.fail()) {
        std::cout << "Data corruption at fpol.\n";
        return is;
    }
    // pressure
    read_vec(g.nw, g.presure);
    if (is.fail()) {
        std::cout << "Data corruption at pressure.\n";
        return is;
    }
    // f*f^{\prime}
    read_vec(g.nw, g.f_f_prime);
    if (is.fail()) {
        std::cout << "Data corruption at ffprime.\n";
        return is;
    }
    // p^{\prime}
    read_vec(g.nw, g.p_prime);
    if (is.fail()) {
        std::cout << "Data corruption at pprime.\n";
        return is;
    }
    // flux
    g.flux.resize({g.nw, g.nh});
    for (unsigned j = 0; j < g.nh; ++j) {
        for (unsigned i = 0; i < g.nw; ++i) { is >> g.flux(i, j); }
    }
    if (is.fail()) {
        std::cout << "Data corruption at psirz.\n";
        return is;
    }
    // safety factor
    read_vec(g.nw, g.safety_factor);
    if (is.fail()) {
        std::cout << "Data corruption at q.\n";
        return is;
    }
    // # boundary point and # limiter point
    is >> g.boundary_num >> g.limiter_num;
    if (is.fail()) {
        std::cout << "Data corruption at boundary number or limiter number.\n";
        return is;
    }
    // boundary points
    for (unsigned i = 0; i < g.boundary_num; ++i) {
        double rr, zz;
        is >> rr >> zz;
        g.boundary.emplace_back(rr, zz);
    }
    if (is.fail()) {
        std::cout << "Data corruption at boundary points.\n";
        return is;
    }

    g.rearrange_boundary();

    // limiter points
    for (unsigned i = 0; i < g.limiter_num; ++i) {
        double rr, zz;
        is >> rr >> zz;
        g.limiter.emplace_back(rr, zz);
    }
    if (is.fail()) {
        std::cout << "Data corruption at limiter points.\n";
        return is;
    }

    if (!is.fail()) { g.__complete = true; }

    return is;
}

void GFileRawData::rearrange_boundary() {
    geometric_poloidal_angles.reserve(boundary_num);

    size_t middle = 0;
    for (size_t i = 0; i < boundary_num; ++i) {
        geometric_poloidal_angles.push_back(
            util::arctan(boundary[i] - magnetic_axis));

        if (i > 0 && std::abs(geometric_poloidal_angles[i] -
                              geometric_poloidal_angles[i - 1]) > M_PI) {
            middle = i;
        }
    }

    std::rotate(boundary.begin(), boundary.begin() + middle, boundary.end());
    std::rotate(geometric_poloidal_angles.begin(),
                geometric_poloidal_angles.begin() + middle,
                geometric_poloidal_angles.end());

    {
        auto last = std::unique(geometric_poloidal_angles.begin(),
                                geometric_poloidal_angles.end());
        geometric_poloidal_angles.erase(last, geometric_poloidal_angles.end());
    }
    {
        auto last = std::unique(boundary.begin(), boundary.end());
        boundary.erase(last, boundary.end());
    }

    // make sure boundary points are sorted counterclockwise
    if (geometric_poloidal_angles[0] > geometric_poloidal_angles[1]) {
        std::reverse(boundary.begin(), boundary.end());
        std::reverse(geometric_poloidal_angles.begin(),
                     geometric_poloidal_angles.end());
    }
}