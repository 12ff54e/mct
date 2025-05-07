#ifndef MEQ_MAGNETIC_EQUILIBRIUM_H
#define MEQ_MAGNETIC_EQUILIBRIUM_H

#include "BSplineInterpolation/src/include/Interpolation.hpp"
#include "GFileRawData.h"

#ifdef MEQ_ZERNIKE_SERIES_
#include "Zernike.h"
#endif

class MagneticEquilibrium {
   public:
    // interpolation order of internal use
    constexpr static std::size_t ORDER = 5;
    constexpr static std::size_t ORDER_OUT = 2;

   protected:
    static constexpr std::size_t FIELD_NUM_2D = 4;
    static constexpr std::size_t FIELD_NUM_1D = 6;

    // The contour mesh grid of data stored here, and the resolution should be
    // higher than the output spec
    struct MagneticEquilibriumRaw_ {
        // - magnetic_field
        // - r
        // - z
        // - jacobian
        std::array<intp::Mesh<double, 2>, FIELD_NUM_2D> data_2d;
        // - safety_factor
        // - poloidal_current
        // - toroidal_current
        // - pressure
        // - minor_radius
        // - toroidal_flux
        std::array<std::vector<double>, FIELD_NUM_1D> data_1d;

        std::array<double, FIELD_NUM_2D> axis_value_2d;
        std::array<double, FIELD_NUM_1D> axis_value_1d;

        double flux_unit;
    };

    struct MagneticEquilibriumIntp_ {
        std::vector<double> psi_sample_for_output;

        // - magnetic_field
        // - r
        // - z
        // - jacobian
#ifdef MEQ_ZERNIKE_SERIES_
        std::array<Zernike::Series<double>, FIELD_NUM_2D>
#else
        std::array<intp::InterpolationFunction<double, 2, ORDER_OUT>,
                   FIELD_NUM_2D>
#endif
            intp_2d;
        // - safety_factor
        // - poloidal_current
        // - toroidal_current
        // - pressure
        // - minor_radius
        // - toroidal_flux
        std::array<intp::InterpolationFunction1D<ORDER_OUT, double>,
                   FIELD_NUM_1D>
            intp_1d;

        template <typename F,
                  std::size_t... indices_2d,
                  std::size_t... indices_1d>
        MagneticEquilibriumIntp_(const MagneticEquilibriumRaw_& spdata_raw,
                                 const MagneticEquilibrium& spdata,
                                 const F& generate_psi_for_output,
                                 std::index_sequence<indices_2d...>,
                                 std::index_sequence<indices_1d...>)
            : psi_sample_for_output{generate_psi_for_output(spdata.psi_delta(),
                                                            spdata.lsp)},
              intp_2d{spdata.create_2d_spline_(spdata_raw.data_2d[indices_2d],
                                               psi_sample_for_output)...},
              intp_1d{spdata.create_1d_spline_(spdata_raw.data_1d[indices_1d],
                                               psi_sample_for_output)...} {}
    };

   public:
    template <typename F>
    MagneticEquilibrium(const GFileRawData& g_file_data,
                        std::size_t radial_grid_num,
                        std::size_t poloidal_grid_num,
                        bool use_si,
                        std::size_t radial_sample,
                        double psi_ratio,
                        const F& generate_psi_for_output)
        : lsp(radial_grid_num),
          lst(poloidal_grid_num),
          use_si_(use_si),
          psi_delta_(psi_ratio *
                     (g_file_data.flux_LCFS - g_file_data.flux_magnetic_axis) /
                     static_cast<double>(lsp - 1)),
          theta_delta_(2. * M_PI / static_cast<double>(lst)),
          spdata_raw_{generate_boozer_coordinate_(g_file_data,
                                                  radial_sample,
                                                  psi_ratio)},
          spdata_intp_{spdata_raw_, *this, generate_psi_for_output,
                       std::make_index_sequence<FIELD_NUM_2D>{},
                       std::make_index_sequence<FIELD_NUM_1D>{}} {}

    const std::size_t lsp, lst;

    const std::vector<double>& psi_for_output() const;
    double psi_delta() const;
    double theta_delta() const;
    const MagneticEquilibriumIntp_& intp_data() const;
    const std::array<double, FIELD_NUM_2D>& axis_value_2d() const;
    const std::array<double, FIELD_NUM_1D>& axis_value_1d() const;

   private:
    const bool use_si_;
    double psi_delta_;
    const double theta_delta_;
    const MagneticEquilibriumRaw_ spdata_raw_;
    const MagneticEquilibriumIntp_ spdata_intp_;

    MagneticEquilibriumRaw_ generate_boozer_coordinate_(const GFileRawData&,
                                                        std::size_t,
                                                        double);
#ifdef MEQ_ZERNIKE_SERIES_
    Zernike::Series<double>
#else
    intp::InterpolationFunction<double, 2, ORDER_OUT>
#endif
    create_2d_spline_(const intp::Mesh<double, 2>&,
                      const std::vector<double>&) const;
    intp::InterpolationFunction1D<ORDER_OUT, double> create_1d_spline_(
        const std::vector<double>&,
        const std::vector<double>&) const;
};

#endif  // MEQ_MAGNETIC_EQUILIBRIUM_H
