#include <ostream>

#include "MagneticEquilibrium.h"

struct Spdata : MagneticEquilibrium<double> {
    Spdata(const GFileRawData<val_type>& g_file_data,
           std::size_t radial_grid_num,
           std::size_t poloidal_grid_num,
           bool use_si = false,
           std::size_t radial_sample = 256,
           val_type psi_ratio = .99);

    void print(std::ostream&) const;

   private:
    static std::vector<val_type> generate_psi_for_output_(val_type,
                                                          std::size_t);
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
