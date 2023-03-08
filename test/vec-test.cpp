#include "../src/include/Vec.hpp"
#include <iomanip>
#include "Assertion.hpp"

int main() {
    Assertion assertion;

    auto v2_origin = Vec<2, double>::zero();
    Vec<2, double> v2_a{1, 2};

    assertion(v2_origin.x() == 0,
              std::string{"Difference between origin x and 0: "} +
                  std::to_string(v2_origin.x()));
    assertion(v2_a.y() == 2, std::string{"Difference between {1,2}.y and 2: "} +
                                 std::to_string(v2_a.y() - 2));

    Vec<3, int> v3_a{2, 2, 2};
    Vec<3, double> v3_b(v3_a);

    assertion(v3_a.x() == 2., "Convert construction failed.");

    auto v3_c = v3_b;
    v3_c.z() = 4;
    v3_b += v3_c;

    assertion(v3_b.x() == 4 && v3_b.y() == 4 && v3_b.z() == 6,
              "Compound assignment by sum failed.");

    auto v2_b = v2_a * 2;
    auto v2_c = 10 * v2_a;
    assertion(
        v2_b.x() == 2 && v2_b.y() == 4 && v2_c.x() == 10 && v2_c.y() == 20,
        "Scalar product failed.");

    auto v2_sum = v2_b + v2_c;
    auto v2_diff = v2_c - v2_b;
    assertion(
        v2_sum == Vec<2, double>{12, 24} && v2_diff == Vec<2, double>{8, 16},
        "Sum or Diff or Comparison failed.");

    assertion(cross(Vec<2>{1, 2}) == Vec<2>{-2, 1}, "Cross on Vec2 failed.");
    assertion(cross(Vec<3>{1, 2, 3}, Vec<3>{1, 2, 3}) == Vec<3>::zero(),
              "Cross product of Vec3 failed.");
    assertion(dot(Vec<3>{1, 2, 3}, Vec<3>{3, 2, 1}) == 10,
              "Dot product failed");

    return assertion.status();
}
