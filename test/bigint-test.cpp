#include <iostream>

#include "../lib/BigInt.h"
#include "Assertion.hpp"

int main() {
    Assertion assertion;

    {
        using i0 = BigInt<Polynomial<std::ratio<499>, std::ratio<1>>>;
        using i1 = BigInt<Polynomial<std::ratio<999>>>;
        assertion.test(!bigint_less_v<i0, i1> && bigint_less_v<i1, i0>,
                       "less test");
    }
    {
        using i0 = BigInt<Polynomial<std::ratio<42>>>;
        using i1 = BigInt<Polynomial<std::ratio<42>>>;
        using i2 = BigInt<Polynomial<std::ratio<69>>>;
        assertion.test(bigint_equal_v<i0, i1> && bigint_equal_v<i1, i0> &&
                           !bigint_equal_v<i0, i2>,
                       "equal test");
        assertion.test(
            bigint_less_equal_v<i0, i1> && bigint_less_equal_v<i0, i2>,
            "less_equal test");
    }
    {
        using i0 = BigInt<Polynomial<std::ratio<499>, std::ratio<1>>>;
        using i1 = BigInt<Polynomial<std::ratio<999>>>;
        assertion.test(
            std::is_same_v<bigint_add<i0, i1>,
                           BigInt<Polynomial<std::ratio<498>, std::ratio<2>>>>,
            "addition test");
    }
    {
        using i0 = BigInt<Polynomial<std::ratio<1>, std::ratio<1>>>;
        using i1 = BigInt<Polynomial<std::ratio<499>>>;
        using i2 =
            BigInt<Polynomial<std::ratio<0>, std::ratio<0>, std::ratio<1>>>;
        assertion.test(
            std::is_same_v<bigint_sub<i0, i1>,
                           BigInt<Polynomial<std::ratio<502>>>> &&
                std::is_same_v<
                    bigint_sub<i0, i2>,
                    BigInt<Polynomial<std::ratio<-999>, std::ratio<-998>>>>,
            "subtraction test");
    }
    {
        using i0 = BigInt<Polynomial<std::ratio<999>, std::ratio<1>>>;
        using i1 = BigInt<Polynomial<std::ratio<999>>>;
        assertion.test(
            std::is_same_v<bigint_mul<i0, i1>,
                           BigInt<Polynomial<std::ratio<1>, std::ratio<997>,
                                             std::ratio<1>>>>,
            "multiplication test");
    }

    // {
    //     using i0 = BigInt<Polynomial<std::ratio<200>, std::ratio<1>>>;
    //     using i1 = BigInt<Polynomial<std::ratio<240>>>;
    //     assertion.test(std::is_same_v<bigint_div<i0, i1>,
    //                                   BigInt<Polynomial<std::ratio<5>>>>,
    //                    "division test");
    // }
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

        //     assertion.test(
        //         std::is_same_v<
        //             bigint_factorial_partial<to_bigint<20>, to_bigint<10>>,
        //             bigint_div<bigint_factorial<to_bigint<20>>,
        //                        bigint_factorial<to_bigint<10>>>>,
        //         "partial factorial test");
    }
    return assertion.status();
}
