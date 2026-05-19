#include <ostream>

#include "MagneticEquilibrium.h"

struct Spdata {
    // NOTE: the only location that decides the value type of everything
    using val_type = double;

    const std::size_t lsp;
    const std::size_t lst;

    Spdata(const GFileRawData<val_type>& g_file_data,
           std::size_t radial_grid_num,
           std::size_t poloidal_grid_num,
           bool use_si = false,
           std::size_t radial_sample = 256,
           val_type psi_ratio = .99);

    void print(std::ostream&) const;

    val_type psi_delta(val_type coef = val_type{1.}) const;
    val_type theta_delta(val_type coef = val_type{1.}) const;

   private:
    constexpr static std::size_t ORDER = 5;
    // output data is coefficient of quadratic spline
    constexpr static std::size_t ORDER_OUT = 2;
    constexpr static std::size_t FIELD_NUM_1D = 6;
    constexpr static std::size_t FIELD_NUM_2D = 4;

    val_type psi_delta_;
    val_type theta_delta_;

    std::array<val_type, FIELD_NUM_1D> axis_value_1d_;
    std::array<val_type, FIELD_NUM_2D> axis_value_2d_;
    // - safety_factor
    // - poloidal_current
    // - toroidal_current
    // - pressure
    // - minor_radius
    // - toroidal_flux
    std::array<intp::InterpolationFunction1D<ORDER_OUT, val_type, val_type>,
               FIELD_NUM_1D>
        intp_1d_{};
    // - magnetic_field
    // - r
    // - z
    // - jacobian
    std::array<intp::InterpolationFunction<val_type, 2, ORDER_OUT, val_type>,
               FIELD_NUM_2D>
        intp_2d_{};

    const std::vector<val_type>& psi_for_output_() const;

    intp::InterpolationFunction<val_type, 2, ORDER_OUT, val_type>
    create_2d_spline_(const intp::Mesh<val_type, 2>& data) const;
};

struct SpdataLiteral {
    using val_type = Spdata::val_type;
    std::string meta_data;

    std::size_t lsp;
    std::size_t lst;
    val_type psi_wall;
    val_type psi_sep;

    // - magnetic_field
    // - r
    // - z
    // - jacobian
    intp::Mesh<val_type, 3> b;
    intp::Mesh<val_type, 3> r;
    intp::Mesh<val_type, 3> z;
    intp::Mesh<val_type, 3> j;

    // - safety_factor
    // - poloidal_current
    // - toroidal_current
    // - pressure
    // - minor_radius
    // - toroidal_flux
    intp::Mesh<val_type, 2> q;
    intp::Mesh<val_type, 2> f;
    intp::Mesh<val_type, 2> i;
    intp::Mesh<val_type, 2> p;
    intp::Mesh<val_type, 2> r_minor;
    intp::Mesh<val_type, 2> psi_t;

    SpdataLiteral(std::size_t, std::size_t);

    static std::tuple<std::size_t, std::size_t, std::size_t> peek_dimension(
        std::istream&);
    static SpdataLiteral construct_from_input(std::istream&);
    void print_to(std::ostream&) const;
    void convert_to_SI(val_type, val_type);

   private:
    void read_from(std::istream&, std::size_t);
};
