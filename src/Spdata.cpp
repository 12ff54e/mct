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

#ifndef MCT_ZERNIKE_SERIES_
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
