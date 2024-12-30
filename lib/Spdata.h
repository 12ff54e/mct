#ifndef MCT_SPDATA_
#define MCT_SPDATA_

#include <ostream>

#include "BSplineInterpolation/src/include/Interpolation.hpp"
#include "GFileRawData.h"

class Spdata {
    static constexpr std::size_t FIELD_NUM_2D = 4;
    static constexpr std::size_t FIELD_NUM_1D = 6;

    // The contour mesh grid of data stored here, and the resolution should be
    // higher than the output spec
    struct SpdataRaw_ {
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

    struct SpdataIntp_ {
        std::vector<double> psi_sample_for_output;

        // - magnetic_field
        // - r
        // - z
        // - jacobian
        std::array<intp::InterpolationFunction<double, 2>, FIELD_NUM_2D>
            intp_2d;
        // - safety_factor
        // - poloidal_current
        // - toroidal_current
        // - pressure
        // - minor_radius
        // - toroidal_flux
        std::array<intp::InterpolationFunction1D<double>, FIELD_NUM_1D> intp_1d;

        template <std::size_t... indices_2d, std::size_t... indices_1d>
        SpdataIntp_(const SpdataRaw_& spdata_raw,
                    const Spdata& spdata,
                    std::index_sequence<indices_2d...>,
                    std::index_sequence<indices_1d...>)
            : psi_sample_for_output{spdata.generate_psi_sample_for_output_(
                  spdata_raw.flux_unit)},
              intp_2d{spdata.create_2d_spline_(spdata_raw.data_2d[indices_2d],
                                               psi_sample_for_output)...},
              intp_1d{spdata.create_1d_spline_(spdata_raw.data_1d[indices_1d],
                                               psi_sample_for_output)...} {}
    };

   public:
    Spdata(const GFileRawData&,
           std::size_t,
           std::size_t,
           bool = false,
           std::size_t = 256,
           double = .99);

    friend std::ostream& operator<<(std::ostream&, const Spdata&);

   private:
    // interpolation order of internal use
    constexpr static std::size_t ORDER_ = 5;
    constexpr static std::size_t ORDER_OUT_ = 2;

    const bool use_si_;
    const std::size_t lsp_, lst_;
    double psi_delta_;
    const double theta_delta_;
    const SpdataRaw_ spdata_raw_;
    const SpdataIntp_ spdata_intp_;

    SpdataRaw_ generate_boozer_coordinate_(const GFileRawData&,
                                           std::size_t,
                                           double);
    std::vector<double> generate_psi_sample_for_output_(double) const;
    intp::InterpolationFunction<double, 2> create_2d_spline_(
        const intp::Mesh<double, 2>&,
        const std::vector<double>&) const;
    intp::InterpolationFunction1D<double> create_1d_spline_(
        const std::vector<double>&,
        const std::vector<double>&) const;
};

#endif
