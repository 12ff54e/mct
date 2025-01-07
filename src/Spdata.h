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
