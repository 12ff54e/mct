#include "../lib/Polynomial.h"

int main() {
    using P1 = Polynomial<std::ratio<1>, std::ratio<0>, std::ratio<1>>;
    using P2 = Polynomial<std::ratio<1>, std::ratio<1>>;

    using quotient = Polynomial<std::ratio<-1>, std::ratio<1>>;
    using remainder = Polynomial<std::ratio<2>>;

    static_assert(std::is_same_v<poly_div<P1, P2>, quotient>);
    static_assert(std::is_same_v<poly_mod<P1, P2>, remainder>);
    return 0;
}
