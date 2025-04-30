#include <iostream>

#include "../lib/BigInt.h"
#include "Assertion.hpp"

int main() {
    Assertion assertion;
    {
        using i0 = BigInt<Polynomial<std::ratio<999>, std::ratio<1>>>;
        using i1 = BigInt<Polynomial<std::ratio<999>>>;
        assertion.test(
            std::is_same_v<bigint_mul<i0, i1>,
                           BigInt<Polynomial<std::ratio<1>, std::ratio<997>,
                                             std::ratio<1>>>>,
            "multiplication test");
    }
    {
        using i0 = BigInt<Polynomial<std::ratio<200>, std::ratio<1>>>;
        using i1 = BigInt<Polynomial<std::ratio<240>>>;
        assertion.test(std::is_same_v<bigint_div<i0, i1>,
                                      BigInt<Polynomial<std::ratio<5>>>>,
                       "division test");
    }
    {
        assertion.test(
            std::is_same_v<to_bigint<1234567>,
                           BigInt<Polynomial<std::ratio<567>, std::ratio<234>,
                                             std::ratio<1>>>>,
            "conversion test");
    }
    {
        assertion.test(
            std::is_same_v<bigint_factorial<to_bigint<10>>,
                           BigInt<Polynomial<std::ratio<800>, std::ratio<628>,
                                             std::ratio<3>>>>,
            "factorial test");
    }
    return assertion.status();
}
