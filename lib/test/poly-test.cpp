#include "Polynomial.h"

int main() {
    using P1 = Polynomial<1, 2>;
    using P2 = Polynomial<3, 4>;

    using product = Polynomial<3, 10, 8>;

    static_assert(std::is_same_v<poly_mul<P1, P2>, product>);
    return 0;
}
