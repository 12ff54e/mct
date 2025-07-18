#include "Spdata.h"
#include <sstream>

#include <iomanip>

Spdata::Spdata(const GFileRawData<val_type>& g_file_data,
               std::size_t radial_grid_num,
               std::size_t poloidal_grid_num,
               bool use_si,
               std::size_t radial_sample,
               val_type psi_ratio)
    : MagneticEquilibrium(g_file_data,
                          radial_grid_num,
                          poloidal_grid_num,
                          use_si,
                          radial_sample,
                          psi_ratio,
                          generate_psi_for_output_) {}

std::vector<Spdata::val_type> Spdata::generate_psi_for_output_(
    val_type psi_delta,
    std::size_t l) {
    std::vector<val_type> psi(l);
    psi[0] = psi_delta;
    for (std::size_t i = 1; i < l - 1; ++i) {
        psi[i] = psi_delta * (static_cast<val_type>(i) + .5);
    }
    psi[l - 1] = static_cast<val_type>(l - 1) * psi_delta;

    return psi;
}

void Spdata::print(std::ostream& output_stream) const {
    // mesh grid info
    output_stream << std::setw(4) << lsp << std::setw(4) << lst << std::setw(4)
                  << 4 << std::setw(4) << 8 << '\n';
    output_stream << std::scientific << std::uppercase << std::setprecision(10);
    output_stream << std::setw(18) << psi_for_output().back() << std::setw(18)
                  << psi_for_output().back() << '\n';

    std::vector<std::ostringstream> output_buffer(lsp);

#ifndef MEQ_ZERNIKE_SERIES_
    auto write_1d_coef = [psi_delta = psi_delta()](
                             auto& os, const auto& f_1d, std::size_t idx,
                             val_type value_on_axis, bool singular) {
        val_type c0, c1, c2;
        if (idx == 0) {
            const val_type v1 = f_1d(psi_delta) - value_on_axis;
            const val_type v2 = f_1d.derivative(std::make_pair(psi_delta, 1)) *
                                (singular ? 2. * std::sqrt(psi_delta) : 1.);
            const val_type d = singular ? std::sqrt(psi_delta) : psi_delta;
            c0 = value_on_axis;
            c1 = 2. * v1 / d - v2;
            c2 = (-v1 + v2 * d) / (d * d);
        } else {
            const val_type psi = static_cast<val_type>(idx) * psi_delta;
            c0 = f_1d(psi);
            c1 = f_1d.derivative(std::make_pair(psi, 1));
            c2 = .5 * f_1d.derivative(std::make_pair(
                          psi + static_cast<val_type>(.5 * psi_delta), 2));
        }
        os << std::scientific << std::uppercase << std::setprecision(10)
           << std::setw(18) << c0 << std::setw(18) << c1 << std::setw(18) << c2
           << '\n';
    };

    auto write_2d_coef = [&](auto& os, auto& f_2d, std::size_t idx,
                             val_type value_on_axis) {
        os << std::scientific << std::uppercase << std::setprecision(10);
        for (size_t i = 0; i < 9; ++i) {
            size_t psi_order = i % 3;
            size_t theta_order = i / 3;
            for (size_t j = 0; j <= lst; ++j) {
                val_type coef =
                    (psi_order == 2 ? .5 : 1.) * (theta_order == 2 ? .5 : 1.);
                const val_type theta =
                    (static_cast<val_type>(j) + (theta_order == 2 ? .5 : 0.)) *
                    theta_delta();
                if (idx == 0) {
                    const val_type v1 =
                        theta_order == 0
                            ? f_2d(psi_delta(), theta) - value_on_axis
                            : f_2d.derivative({psi_delta(), theta},
                                              {0, theta_order});
                    const val_type v2 =
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
                    const val_type psi = (static_cast<val_type>(idx) +
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

#ifdef MCT_MULTITHREAD
    auto& thread_pool = intp::DedicatedThreadPool<void>::get_instance();
    std::vector<std::future<void>> tasks;
#endif
    constexpr std::size_t task_size = 2;
    for (std::size_t ri = 0; ri < (lsp + task_size - 1) / task_size; ++ri) {
        const auto start = ri * task_size;
        const auto finish = start + task_size > lsp ? lsp : start + task_size;
        auto write_to_buffer = [&, start, finish]() {
            for (std::size_t i = start; i < finish; ++i) {
                for (std::size_t i_2d = 0; i_2d < FIELD_NUM_2D; ++i_2d) {
                    write_2d_coef(output_buffer[i], intp_data().intp_2d[i_2d],
                                  i, axis_value_2d()[i_2d]);
                }
                for (std::size_t i_1d = 0; i_1d < FIELD_NUM_1D; ++i_1d) {
                    write_1d_coef(output_buffer[i], intp_data().intp_1d[i_1d],
                                  i, axis_value_1d()[i_1d], i_1d == 4);
                }
            }
        };
#ifdef MCT_MULTITHREAD
        tasks.push_back(thread_pool.queue_task(write_to_buffer));
#else
        write_to_buffer();
#endif
    }

#ifdef MCT_MULTITHREAD
    for (auto& res : tasks) { res.get(); }
#endif
    for (auto& oss : output_buffer) { output_stream << oss.str(); }
#endif  // #ifndef MEQ_ZERNIKE_SERIES_

    // TODO: ripple related
    output_stream << std::setw(4) << 0 << std::setw(4) << 0 << '\n';
    output_stream << std::setw(18) << axis_value_2d()[1] << std::setw(18) << 0.
                  << std::setw(18) << 0. << '\n';
    output_stream << std::setw(18) << 0. << std::setw(18) << 0. << '\n';
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
    val_type dump_float;
    is >> lsp >> lst >> dump_int >> dump_int;
    is >> dump_float >> dump_float;

    for (std::size_t t = 0; t < lst; ++t) { is >> dump_float; }
    std::string line;

    // actually number of pts on flux surface, could be lst or lst+1
    std::size_t theta_pts = lst;
    std::getline(is, line);
    {
        std::istringstream iss(line);
        if ((iss >> dump_float)) {
            theta_pts++;
        } else {
            std::getline(is, line);
            iss.clear();
            iss.str(line);
            if (!(iss >> dump_float >> dump_float)) { theta_pts++; }
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

    val_type dump_float;
    auto read_2d_flux = [&](auto& field, auto psi) {
        for (std::size_t c = 0; c < 9; ++c) {
            for (std::size_t t = 0; t < lst; ++t) { is >> field(psi, c, t); }
            if (theta_pts != lst) { is >> dump_float; }
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

void SpdataLiteral::convert_to_SI(val_type b0, val_type r0) {
    const auto psi_unit_inv = 1. / (r0 * r0 * b0);
    const auto psi_unit_2_inv = psi_unit_inv * psi_unit_inv;

    auto convert_2d = [&](auto& field, val_type unit) {
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

    auto convert_1d = [&](auto& field, val_type unit, bool singular = false) {
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
