#include <iostream>

#include "include/GFileRawData.hpp"

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
    for (unsigned i = 0; i < g.nw; ++i) {
        for (unsigned j = 0; j < g.nh; ++j) { is >> g.flux(i, j); }
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
    is >> g.bd_num >> g.limiter_num;
    if (is.fail()) {
        std::cout << "Data corruption at boundary number or limiter number.\n";
        return is;
    }
    // boundary points
    std::vector<double> r_dum;
    std::vector<double> z_dum;
    read_vec(g.bd_num, r_dum);
    read_vec(g.bd_num, z_dum);
    for (unsigned i = 0; i < g.bd_num; ++i) {
        g.boundary.emplace_back(r_dum[i], z_dum[i]);
    }
    if (is.fail()) {
        std::cout << "Data corruption at boundary points.\n";
        return is;
    }
    // calculate anchor points
    unsigned r_min{}, r_max{}, z_min{}, z_max{};
    for (unsigned i = 0; i < g.bd_num; ++i) {
        r_min = r_dum[i] < r_dum[r_min] ? i : r_min;
        r_max = r_dum[i] > r_dum[r_max] ? i : r_max;
        z_min = z_dum[i] < z_dum[z_min] ? i : z_min;
        z_max = z_dum[i] > z_dum[z_max] ? i : z_max;
    }
    g.right_anchor = g.boundary[r_max];
    g.left_anchor = g.boundary[r_min];
    g.top_anchor = g.boundary[z_max];
    g.bottom_anchor = g.boundary[z_min];

    // limiter points
    r_dum.clear();
    z_dum.clear();
    read_vec(g.limiter_num, r_dum);
    read_vec(g.limiter_num, z_dum);
    for (unsigned i = 0; i < g.bd_num; ++i) {
        g.limiter.emplace_back(r_dum[i], z_dum[i]);
    }
    if (is.fail()) {
        std::cout << "Data corruption at limiter points.\n";
        return is;
    }

    if (!is.fail()) { g.__complete = true; }

    return is;
}