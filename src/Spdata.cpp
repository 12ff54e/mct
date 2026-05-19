#include <array>
#include <iomanip>
#include <sstream>

#include "Spdata.h"

Spdata::Spdata(const GFileRawData<val_type>& gfile_data,
               std::size_t radial_grid_num,
               std::size_t poloidal_grid_num,
               bool use_si,
               std::size_t radial_sample,
               val_type psi_ratio)
    : lsp(radial_grid_num), lst(poloidal_grid_num) {
    MagneticEquilibrium equilibrium(
        gfile_data, radial_grid_num, poloidal_grid_num, radial_sample,
        psi_ratio, use_si,
        [](auto raw_values_1d) {
            [[maybe_unused]] auto [psi, f, q, p, b2j] = raw_values_1d;

            return std::array{psi, q, f, p};
        },
        [](auto raw_values_1d, auto raw_values_2d) {
            [[maybe_unused]] auto [psi, f, q, p, b2j] = raw_values_1d;
            [[maybe_unused]] auto [dpdr, dpdz, r, z] = raw_values_2d;

            auto b2 = (f * f + dpdr * dpdr + dpdz * dpdz) / (r * r);
            auto j = b2j / b2;
            // B_p^2*J is used for calculating toroidal current
            auto bp2j = (dpdr * dpdr + dpdz * dpdz) * j / (r * r);

            return std::array{std::sqrt(b2), r, z, j, bp2j};
        });

    psi_delta_ = equilibrium.psi_delta();
    theta_delta_ = equilibrium.theta_delta();
    axis_value_1d_[0] = equilibrium.axis_value_1d(1);
    axis_value_1d_[1] = equilibrium.axis_value_1d(2);
    axis_value_1d_[2] = 0.;
    axis_value_1d_[3] = equilibrium.axis_value_1d(3);
    axis_value_1d_[4] = 0.;
    axis_value_1d_[5] = 0.;
    for (std::size_t i = 0; i < FIELD_NUM_2D; ++i) {
        axis_value_2d_[i] = equilibrium.axis_value_2d(i);
    }

    const std::size_t eq_np = equilibrium.data_2d(0).dim_size(0);
    const std::size_t eq_nt = equilibrium.data_2d(0).dim_size(1);
    intp::InterpolationFunctionTemplate1D<ORDER, val_type, val_type>
        poloidal_template(
            std::make_pair(theta_delta(.5), 2. * M_PI + theta_delta(.5)), eq_nt,
            true);
    // calculate toroidal current
    const auto& bp2j = equilibrium.data_2d(4);
    std::vector<val_type> toroidal_current;
    toroidal_current.reserve(eq_np);
    for (std::size_t i = 0; i < eq_np; ++i) {
        auto bp2j_func = poloidal_template.interpolate(
            std::make_pair(bp2j.begin(1, {i, 0}), bp2j.end(1, {i, 0})));
        val_type integration{};
        for (std::size_t j = 0; j < eq_nt; ++j) {
            integration += util::integrate_coarse(bp2j_func, theta_delta(j),
                                                  theta_delta(j + 1));
        }
        toroidal_current.push_back(integration / (2. * M_PI));
    }

    // create 2d interpolation functions for output
    for (std::size_t i = 0; i < FIELD_NUM_2D; ++i) {
        intp_2d_[i] = create_2d_spline_(equilibrium.data_2d(i));
    }

    // minor_r is defined as R-R_0 at outer mid-plane and always normalized by
    // R0
    std::vector<val_type> minor_r;
    minor_r.reserve(eq_np);
    auto r0 = equilibrium.axis_value_2d(1);
    for (std::size_t i = 0; i < eq_np; ++i) {
        minor_r.push_back(intp_2d_[1](equilibrium.data_1d(0)[i], 0.) / r0 - 1.);
    }

    // create 1d interpolation functions for output
    intp::InterpolationFunctionTemplate1D<ORDER, val_type, val_type>
        radial_template(
            std::make_pair(psi_for_output_().front(), psi_for_output_().back()),
            eq_np);
    intp::InterpolationFunctionTemplate1D<ORDER_OUT, val_type, val_type>
        radial_template_for_output(intp::util::get_range(psi_for_output_()),
                                   psi_for_output_().size());

    std::vector<val_type> data_resampled;
    data_resampled.reserve(psi_for_output_().size());

    std::array<decltype(intp::util::get_range(equilibrium.data_1d(0))), 6>
        data_1d_range;
    data_1d_range[0] = intp::util::get_range(equilibrium.data_1d(1));
    data_1d_range[1] = intp::util::get_range(equilibrium.data_1d(2));
    data_1d_range[2] = intp::util::get_range(std::as_const(toroidal_current));
    data_1d_range[3] = intp::util::get_range(equilibrium.data_1d(3));
    data_1d_range[4] = intp::util::get_range(std::as_const(minor_r));
    // use as_const since equilibrium.data_1d() returns const auto&

    std::vector<val_type> psi_t;
    psi_t.reserve(eq_np);
    // trapezoidal 1st segment
    psi_t.push_back(psi_delta(
        .5 * (equilibrium.axis_value_1d(1) + equilibrium.data_1d(1).front())));

    for (std::size_t i = 0; i < FIELD_NUM_1D; ++i) {
        data_resampled.clear();
        auto data_intp = radial_template.interpolate(data_1d_range[i]);
        for (const auto& psi : psi_for_output_()) {
            data_resampled.push_back(data_intp(psi));
        }
        intp_1d_[i] = radial_template_for_output.interpolate(
            intp::util::get_range(data_resampled));
        // i==0 is safety factor, and from its integration we have psi_t
        if (i == 0) {
            for (std::size_t j = 0; j < eq_np - 1; ++j) {
                psi_t.push_back(
                    psi_t.back() +
                    util::integrate_coarse(data_intp, equilibrium.data_1d(0)[j],
                                           equilibrium.data_1d(0)[j + 1]));
            }
            data_1d_range[5] = intp::util::get_range(std::as_const(psi_t));
        }
    }
}

void Spdata::print(std::ostream& output_stream) const {
    // mesh grid info
    output_stream << std::setw(4) << lsp << std::setw(4) << lst << std::setw(4)
                  << 4 << std::setw(4) << 8 << '\n';
    output_stream << std::scientific << std::uppercase << std::setprecision(10);
    output_stream << std::setw(18) << psi_for_output_().back() << std::setw(18)
                  << psi_for_output_().back() << '\n';

    std::vector<std::ostringstream> output_buffer(lsp);

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
                    write_2d_coef(output_buffer[i], intp_2d_[i_2d], i,
                                  axis_value_2d_[i_2d]);
                }
                for (std::size_t i_1d = 0; i_1d < FIELD_NUM_1D; ++i_1d) {
                    write_1d_coef(output_buffer[i], intp_1d_[i_1d], i,
                                  axis_value_1d_[i_1d], i_1d == 4);
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

    // TODO: ripple related
    output_stream << std::setw(4) << 0 << std::setw(4) << 0 << '\n';
    output_stream << std::setw(18) << axis_value_2d_[1] << std::setw(18) << 0.
                  << std::setw(18) << 0. << '\n';
    output_stream << std::setw(18) << 0. << std::setw(18) << 0. << '\n';
}

Spdata::val_type Spdata::psi_delta(val_type coef) const {
    return psi_delta_ * coef;
}
Spdata::val_type Spdata::theta_delta(val_type coef) const {
    return theta_delta_ * coef;
}

const std::vector<Spdata::val_type>& Spdata::psi_for_output_() const {
    static std::vector<val_type> psi{([&]() {
        std::vector<val_type> p(lsp);
        p[0] = psi_delta();
        for (std::size_t i = 1; i < lsp - 1; ++i) {
            p[i] = psi_delta(static_cast<val_type>(i) + .5);
        }
        p[lsp - 1] = psi_delta(static_cast<val_type>(lsp - 1));
        return p;
    })()};
    return psi;
}

intp::InterpolationFunction<Spdata::val_type,
                            2,
                            Spdata::ORDER_OUT,
                            Spdata::val_type>
Spdata::create_2d_spline_(const intp::Mesh<val_type, 2>& data) const {
    // interpolate the even-spaced data
    intp::InterpolationFunction<val_type, 2, ORDER_OUT, val_type> data_intp(
        {false, true}, data,
        std::make_pair(psi_for_output_().front(), psi_for_output_().back()),
        std::make_pair(theta_delta(.5), 2. * M_PI + theta_delta(.5)));

    // resample on the interpolated function
    intp::Mesh<val_type, 2> data_resampled(lsp, lst + 1);
    for (std::size_t i = 0; i < lsp; ++i) {
        for (std::size_t j = 0; j <= lst; ++j) {
            data_resampled(i, j) =
                data_intp(psi_for_output_()[i],
                          theta_delta(static_cast<val_type>(j) + .5));
        }
    }
    // construct the interpolation function for output
    return {{false, true},
            std::move(data_resampled),
            intp::util::get_range(psi_for_output_()),
            std::make_pair(theta_delta(.5), 2. * M_PI + theta_delta(.5))};
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
