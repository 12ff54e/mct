#include <iostream>

#include "Assertion.hpp"
#include "BigInt.h"

int main() {
    Assertion assertion;

    {
        using i0 = BigInt<Polynomial<499, 1>>;
        using i1 = BigInt<Polynomial<999>>;
        assertion.test(!bigint_less_v<i0, i1> && bigint_less_v<i1, i0>,
                       "less test");
    }
    {
        using i0 = BigInt<Polynomial<42>>;
        using i1 = BigInt<Polynomial<42>>;
        using i2 = BigInt<Polynomial<69>>;
        assertion.test(bigint_equal_v<i0, i1> && bigint_equal_v<i1, i0> &&
                           !bigint_equal_v<i0, i2>,
                       "equal test");
        assertion.test(
            bigint_less_equal_v<i0, i1> && bigint_less_equal_v<i0, i2>,
            "less_equal test");
    }
    {
        using i0 = BigInt<Polynomial<499, 1>>;
        using i1 = BigInt<Polynomial<999>>;
        assertion.test(
            std::is_same_v<bigint_add<i0, i1>, BigInt<Polynomial<498, 2>>>,
            "addition test");
    }
    {
        using i0 = BigInt<Polynomial<1, 1>>;
        using i1 = BigInt<Polynomial<499>>;
        using i2 = BigInt<Polynomial<0, 0, 1>>;
        assertion.test(
            std::is_same_v<bigint_sub<i0, i1>, BigInt<Polynomial<502>>> &&
                std::is_same_v<bigint_sub<i0, i2>,
                               BigInt<Polynomial<-999, -998>>>,
            "subtraction test");
    }
    {
        using i0 = BigInt<Polynomial<999, 1>>;
        using i1 = BigInt<Polynomial<999>>;
        assertion.test(
            std::is_same_v<bigint_mul<i0, i1>, BigInt<Polynomial<1, 997, 1>>>,
            "multiplication test");
    }

    {
        using i0 = BigInt<Polynomial<200, 1>>;
        using i1 = BigInt<Polynomial<240>>;
        assertion.test(
            std::is_same_v<bigint_div<i0, i1>, BigInt<Polynomial<5>>>,
            "division test");
    }
    {
        assertion.test(
            std::is_same_v<to_bigint<1234567>, BigInt<Polynomial<567, 234, 1>>>,
            "int -> BigInt test");
        assertion.test(
            from_bigint_v<BigInt<Polynomial<567, 234, 1>>> == 1234567,
            "BigInt -> int test");
    }
    {
        assertion.test(std::is_same_v<bigint_factorial<to_bigint<10>>,
                                      BigInt<Polynomial<800, 628, 3>>>,
                       "factorial test");

        assertion.test(
            std::is_same_v<
                bigint_factorial_partial<to_bigint<20>, to_bigint<10>>,
                bigint_div<bigint_factorial<to_bigint<20>>,
                           bigint_factorial<to_bigint<10>>>>,
            "partial factorial test");
    }
    return assertion.status();
}
