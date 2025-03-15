#include <ostream>

#include "lib/MagneticEquilibrium.h"

struct Spdata : MagneticEquilibrium {
    Spdata(const GFileRawData& g_file_data,
           std::size_t radial_grid_num,
           std::size_t poloidal_grid_num,
           bool use_si = false,
           std::size_t radial_sample = 256,
           double psi_ratio = .99);

    void print(std::ostream&) const;

   private:
    static std::vector<double> generate_psi_for_output_(double, std::size_t);
};

struct SpdataLiteral {
    std::string meta_data;

    std::size_t lsp;
    std::size_t lst;
    double psi_wall;
    double psi_sep;

    // - magnetic_field
    // - r
    // - z
    // - jacobian
    intp::Mesh<double, 3> b;
    intp::Mesh<double, 3> r;
    intp::Mesh<double, 3> z;
    intp::Mesh<double, 3> j;

    // - safety_factor
    // - poloidal_current
    // - toroidal_current
    // - pressure
    // - minor_radius
    // - toroidal_flux
    intp::Mesh<double, 2> q;
    intp::Mesh<double, 2> f;
    intp::Mesh<double, 2> i;
    intp::Mesh<double, 2> p;
    intp::Mesh<double, 2> r_minor;
    intp::Mesh<double, 2> psi_t;

    SpdataLiteral(std::size_t, std::size_t);

    static std::tuple<std::size_t, std::size_t, std::size_t> peek_dimension(
        std::istream&);
    static SpdataLiteral construct_from_input(std::istream&);
    void print_to(std::ostream&) const;
    void convert_to_SI(double, double);

   private:
    void read_from(std::istream&, std::size_t);
};
