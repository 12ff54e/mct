#include <fstream>
#include <string>
#include <vector>

#include "Vec.hpp"
#include "lib/BSplineInterpolation/intp"

struct GFileRawData {
    std::string code_name;
    std::string date;
    std::string shot_id;
    std::string time;

    // horizontal sample number
    unsigned int nw;
    // vertical sample number
    unsigned int nh;
    // horizontal and vertical dimension
    Vec<2, double> dim;
    // center position of horizontal range
    double r_center;
    // left position of horizontal range
    double r_left;
    // center position of vertical range
    double z_mid;
    // magnetic axis coordinate
    Vec<2, double> magnetic_axis;
    // magnetic flux at magnetic axis
    double flux_magnetic_axis;
    // magnetic flux at Last Closed Flux Surface
    double flux_LCFS;
    // mangetic field strength at center position
    double b_center;
    // current
    double current;
    // magnetic flux at seperatrix
    double flux_sep;
    // seperatrix coordinate
    Vec<2, double> sep;

    // poloidal current density profile
    std::vector<double> f_pol;
    // presure profile
    std::vector<double> presure;
    // f*f^{\prime}
    std::vector<double> f_f_prime;
    // p^{\prime}
    std::vector<double> p_prime;

    // magnetic flux
    intp::Mesh<double, 2> flux;

    // safety factor profile
    std::vector<double> safety_factor;

    // boundary point number
    unsigned int bd_num;
    // limiter pointer number
    unsigned int limiter_num;

    // boundary point coordinates
    std::vector<Vec<2, double>> boundary;
    // limiter point coordinates
    std::vector<Vec<2, double>> limiter;

    // Following data are not read directly from gfile.

    // Anchor points of boundary sample points, might be slightly closer to
    // center than that real boundary.

    Vec<2, double> right_anchor;
    Vec<2, double> left_anchor;
    Vec<2, double> top_anchor;
    Vec<2, double> bottom_anchor;

    // constructors

    GFileRawData() : flux{1, 1}, __complete(false){};
    GFileRawData(const GFileRawData&) = delete;
    GFileRawData(GFileRawData&&) = default;

    // properties

    bool is_complete() const noexcept;

    friend std::ifstream& operator>>(std::ifstream&, GFileRawData&);

   private:
    // the completeness of raw data
    bool __complete;
};
