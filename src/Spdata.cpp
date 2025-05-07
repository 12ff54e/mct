#include "Spdata.h"

#include <iomanip>

Spdata::Spdata(const GFileRawData& g_file_data,
               std::size_t radial_grid_num,
               std::size_t poloidal_grid_num,
               bool use_si,
               std::size_t radial_sample,
               double psi_ratio)
    : MagneticEquilibrium(g_file_data,
                          radial_grid_num,
                          poloidal_grid_num,
                          use_si,
                          radial_sample,
                          psi_ratio,
                          generate_psi_for_output_) {}

std::vector<double> Spdata::generate_psi_for_output_(double psi_delta,
                                                     std::size_t l) {
    std::vector<double> psi(l);
    psi[0] = psi_delta;
    for (std::size_t i = 1; i < l - 1; ++i) {
        psi[i] = psi_delta * (static_cast<double>(i) + .5);
    }
    psi[l - 1] = static_cast<double>(l - 1) * psi_delta;

    return psi;
}

void Spdata::print(std::ostream& os) const {
    // mesh grid info
    os << std::setw(4) << lsp << std::setw(4) << lst << std::setw(4) << 4
       << std::setw(4) << 8 << '\n';
    os << std::scientific << std::uppercase << std::setprecision(10);
    os << std::setw(18) << psi_for_output().back() << std::setw(18)
       << psi_for_output().back() << '\n';

#ifndef MEQ_ZERNIKE_SERIES_
    auto write_1d_coef = [&](const auto& f_1d, std::size_t idx,
                             double value_on_axis, bool singular) {
        double c0, c1, c2;
        if (idx == 0) {
            const double v1 = f_1d(psi_delta()) - value_on_axis;
            const double v2 = f_1d.derivative(std::make_pair(psi_delta(), 1)) *
                              (singular ? 2. * std::sqrt(psi_delta()) : 1.);
            const double d = singular ? std::sqrt(psi_delta()) : psi_delta();
            c0 = value_on_axis;
            c1 = 2. * v1 / d - v2;
            c2 = (-v1 + v2 * d) / (d * d);
        } else {
            const double psi = static_cast<double>(idx) * psi_delta();
            c0 = f_1d(psi);
            c1 = f_1d.derivative(std::make_pair(psi, 1));
            c2 =
                .5 * f_1d.derivative(std::make_pair(psi + .5 * psi_delta(), 2));
        }
        os << std::setw(18) << c0 << std::setw(18) << c1 << std::setw(18) << c2
           << '\n';
    };

    auto write_2d_coef = [&](auto& f_2d, std::size_t idx,
                             double value_on_axis) {
        for (size_t i = 0; i < 9; ++i) {
            size_t psi_order = i % 3;
            size_t theta_order = i / 3;
            for (size_t j = 0; j <= lst; ++j) {
                double coef =
                    (psi_order == 2 ? .5 : 1.) * (theta_order == 2 ? .5 : 1.);
                const double theta =
                    (static_cast<double>(j) + (theta_order == 2 ? .5 : 0.)) *
                    theta_delta();
                if (idx == 0) {
                    const double v1 =
                        theta_order == 0
                            ? f_2d(psi_delta(), theta) - value_on_axis
                            : f_2d.derivative({psi_delta(), theta},
                                              {0, theta_order});
                    const double v2 =
                        f_2d.derivative({psi_delta(), theta}, {1, theta_order});

                    coef = psi_order == 0
                               ? theta_order == 0 ? value_on_axis : 0.
                               : 2. * coef *
                                     (psi_order == 1
                                          ? (v1 - v2 * psi_delta()) /
                                                std::sqrt(psi_delta())
                                          : -(v1 - 2. * v2 * psi_delta()) /
                                                psi_delta());
                } else {
                    const double psi = (static_cast<double>(idx) +
                                        (psi_order == 2 ? .5 : 0.)) *
                                       psi_delta();
                    coef = (i == 0 ? f_2d(psi, theta)
                                   : coef * f_2d.derivative(
                                                {psi, theta},
                                                {psi_order, theta_order}));
                }

                os << std::setw(18) << coef;
                if (j % 4 == 3) { os << '\n'; }
            }
            if ((lst + 1) % 4 != 0) { os << '\n'; }
        }
    };

    for (std::size_t ri = 0; ri < lsp; ++ri) {
        for (std::size_t i_2d = 0; i_2d < FIELD_NUM_2D; ++i_2d) {
            write_2d_coef(intp_data().intp_2d[i_2d], ri, axis_value_2d()[i_2d]);
        }
        for (std::size_t i_1d = 0; i_1d < FIELD_NUM_1D; ++i_1d) {
            write_1d_coef(intp_data().intp_1d[i_1d], ri, axis_value_1d()[i_1d],
                          i_1d == 4);
        }
    }
#endif

    // TODO: ripple related
    os << std::setw(4) << 0 << std::setw(4) << 0 << '\n';
    os << std::setw(18) << axis_value_2d()[1] << std::setw(18) << 0.
       << std::setw(18) << 0. << '\n';
    os << std::setw(18) << 0. << std::setw(18) << 0. << '\n';
}

SpdataLiteral::SpdataLiteral(std::size_t psi_num, std::size_t theta_num)
    : lsp(psi_num),
      lst(theta_num),
      b(lsp, 9, lst),
      r(lsp, 9, lst),
      z(lsp, 9, lst),
      j(lsp, 9, lst),
      q(lsp, 3),
      f(lsp, 3),
      i(lsp, 3),
      p(lsp, 3),
      r_minor(lsp, 3),
      psi_t(lsp, 3) {}

std::tuple<std::size_t, std::size_t, std::size_t> SpdataLiteral::peek_dimension(
    std::istream& is) {
    auto max_val = std::numeric_limits<std::streamsize>::max();
    is.ignore(max_val, '\n');

    std::size_t lsp, lst;
    int dump_int;
    double dump_double;
    is >> lsp >> lst >> dump_int >> dump_int;
    is >> dump_double >> dump_double;

    for (std::size_t t = 0; t < lst; ++t) { is >> dump_double; }
    std::string line;

    // actually number of pts on flux surface, could be lst or lst+1
    std::size_t theta_pts = lst;
    std::getline(is, line);
    {
        std::istringstream iss(line);
        if ((iss >> dump_double)) {
            theta_pts++;
        } else {
            std::getline(is, line);
            iss.clear();
            iss.str(line);
            if (!(iss >> dump_double >> dump_double)) { theta_pts++; }
        }
    }

    return {lsp, lst, theta_pts};
}

SpdataLiteral SpdataLiteral::construct_from_input(std::istream& is) {
    auto [lsp, lst, theta_pts] = peek_dimension(is);
    SpdataLiteral spdata(lsp, lst);

    is.seekg(0);
    spdata.read_from(is, theta_pts);

    return spdata;
}

void SpdataLiteral::read_from(std::istream& is, std::size_t theta_pts) {
    std::getline(is, meta_data);

    const auto max_val = std::numeric_limits<std::streamsize>::max();
    is.ignore(max_val, '\n');
    is >> psi_wall >> psi_sep;

    double dump_double;
    auto read_2d_flux = [&](auto& field, auto psi) {
        for (std::size_t c = 0; c < 9; ++c) {
            for (std::size_t t = 0; t < lst; ++t) { is >> field(psi, c, t); }
            if (theta_pts != lst) { is >> dump_double; }
        }
    };

    for (std::size_t psi = 0; psi < lsp; ++psi) {
        read_2d_flux(b, psi);
        read_2d_flux(r, psi);
        read_2d_flux(z, psi);
        read_2d_flux(j, psi);

        is >> q(psi, 0) >> q(psi, 1) >> q(psi, 2);
        is >> f(psi, 0) >> f(psi, 1) >> f(psi, 2);
        is >> i(psi, 0) >> i(psi, 1) >> i(psi, 2);
        is >> p(psi, 0) >> p(psi, 1) >> p(psi, 2);
        is >> r_minor(psi, 0) >> r_minor(psi, 1) >> r_minor(psi, 2);
        is >> psi_t(psi, 0) >> psi_t(psi, 1) >> psi_t(psi, 2);
    }
}

void SpdataLiteral::print_to(std::ostream& os) const {
    os << meta_data << '\n';
    os << std::setw(4) << lsp << std::setw(4) << lst << std::setw(4) << 4
       << std::setw(4) << 8 << '\n';
    os << std::scientific << std::uppercase << std::setprecision(10);
    os << std::setw(18) << psi_wall << std::setw(18) << psi_sep << '\n';

    auto print_2d_flux = [&](const auto& field, std::size_t psi) {
        for (std::size_t c = 0; c < 9; ++c) {
            for (std::size_t t = 0; t < lst + 1; ++t) {
                os << std::setw(18) << field(psi, c, t % lst);
                if (t % 4 == 3) { os << '\n'; }
            }
            if (lst % 4 != 3) { os << '\n'; }
        }
    };

    for (std::size_t psi = 0; psi < lsp; ++psi) {
        print_2d_flux(b, psi);
        print_2d_flux(r, psi);
        print_2d_flux(z, psi);
        print_2d_flux(j, psi);

        os << std::setw(18) << q(psi, 0) << std::setw(18) << q(psi, 1)
           << std::setw(18) << q(psi, 2) << '\n';
        os << std::setw(18) << f(psi, 0) << std::setw(18) << f(psi, 1)
           << std::setw(18) << f(psi, 2) << '\n';
        os << std::setw(18) << i(psi, 0) << std::setw(18) << i(psi, 1)
           << std::setw(18) << i(psi, 2) << '\n';
        os << std::setw(18) << p(psi, 0) << std::setw(18) << p(psi, 1)
           << std::setw(18) << p(psi, 2) << '\n';
        os << std::setw(18) << r_minor(psi, 0) << std::setw(18)
           << r_minor(psi, 1) << std::setw(18) << r_minor(psi, 2) << '\n';
        os << std::setw(18) << psi_t(psi, 0) << std::setw(18) << psi_t(psi, 1)
           << std::setw(18) << psi_t(psi, 2) << '\n';
    }
    os << std::setw(4) << 0 << std::setw(4) << 0 << '\n';
    os << std::setw(18) << 1. << std::setw(18) << 0. << std::setw(18) << 0.
       << '\n';
    os << std::setw(18) << 0. << std::setw(18) << 0. << '\n';
}

void SpdataLiteral::convert_to_SI(double b0, double r0) {
    const auto psi_unit_inv = 1. / (r0 * r0 * b0);
    const auto psi_unit_2_inv = psi_unit_inv * psi_unit_inv;

    auto convert_2d = [&](auto& field, double unit) {
        for (std::size_t psi = 0; psi < lsp; ++psi) {
            for (std::size_t c = 0; c < 9; ++c) {
                for (std::size_t t = 0; t < lst; ++t) {
                    field(psi, c, t) *=
                        unit *
                        (c % 3 == 0 ? 1.
                         : c % 3 == 1
                             ? (psi == 0 ? std::sqrt(psi_unit_inv)
                                         : psi_unit_inv)
                             : (psi == 0 ? psi_unit_inv : psi_unit_2_inv));
                }
            }
        }
    };

    auto convert_1d = [&](auto& field, double unit, bool singular = false) {
        for (std::size_t psi = 0; psi < lsp; ++psi) {
            for (std::size_t c = 0; c < 3; ++c) {
                field(psi, c) *=
                    unit *
                    (c == 0   ? 1.
                     : c == 1 ? (psi == 0 && singular ? std::sqrt(psi_unit_inv)
                                                      : psi_unit_inv)
                              : (psi == 0 && singular ? psi_unit_inv
                                                      : psi_unit_2_inv));
            }
        }
    };

    psi_wall /= psi_unit_inv;
    psi_sep /= psi_unit_inv;

    convert_2d(b, b0);
    convert_2d(r, r0);
    convert_2d(z, r0);
    convert_2d(j, r0 / b0);

    convert_1d(q, 1.);
    convert_1d(f, r0 * b0);
    convert_1d(i, r0 * b0);
    convert_1d(p, 1.);
    convert_1d(r_minor, 1., true);
    convert_1d(psi_t, 1. / psi_unit_inv);
}
