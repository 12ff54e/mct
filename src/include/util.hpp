#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace util {

/**
 * @brief Sign function
 *
 * @return sign of val
 */
template <typename Tf>
int sgn(const Tf& val) {
    return (Tf(0) < val) - (val < Tf(0));
}

template <typename T>
T abs(T x) {
    return std::abs(x);
}

namespace detail {

template <typename Func, typename Tx, typename Tf>
void bracket(Func& func, Tx& a, Tx& b, Tx c, Tf& fa, Tf& fb, Tf& d, Tf& fd) {
    Tf tol = std::numeric_limits<Tf>::epsilon() * 2;

    if ((b - a) < 2 * tol * a) {
        c = a + (b - a) / 2;
    } else if (c - a <= util::abs(a) * tol) {
        c = a + util::abs(a) * tol;
    } else if (b - c <= util::abs(b) * tol) {
        c = b - util::abs(b) * tol;
    }

    Tf fc = func(c);
    if (fc == 0) {
        a = c;
        fa = 0;
        d = 0;
        fd = 0;
        return;
    }

    if (sgn(fa) * sgn(fc) < 0) {
        d = b;
        fd = fb;
        b = c;
        fb = fc;
    } else {
        d = a;
        fd = fa;
        a = c;
        fa = fc;
    }
}

template <typename T>
T safe_div(T num, T denom, T r) {
    return util::abs(denom) < 1 &&
                   util::abs(denom * std::numeric_limits<T>::max()) <=
                       util::abs(num)
               ? r
               : num / denom;
}

template <typename Tx, typename Tf>
Tx secant_interpolate(const Tx& a, const Tx& b, const Tf& fa, const Tf& fb) {
    Tf tol = std::numeric_limits<Tf>::epsilon() * 5;
    Tx c = a - (fa / (fb - fa)) * (b - a);

    return c - a <= util::abs(a) * tol || b - c <= util::abs(b) * tol
               ? (a + b) / 2
               : c;
}

template <typename T>
bool outside(const T& c, const T& a, const T& b) {
    return c <= a || c >= b;
}

template <typename Tx, typename Tf>
Tx quadratic_interpolation(const Tx& a,
                           const Tx& b,
                           const Tx& d,
                           const Tf& fa,
                           const Tf& fb,
                           const Tf& fd,
                           unsigned count) {
    Tf B = safe_div(Tf(fb - fa), Tf(b - a), std::numeric_limits<Tf>::max());
    Tf A = safe_div(Tf(fd - fb), Tf(d - b), std::numeric_limits<Tf>::max());
    A = safe_div(Tf(A - B), Tf(d - a), Tf(0));

    if (A == 0) { return secant_interpolate(a, b, fa, fb); }

    Tx c = sgn(A) * sgn(fa) > 0 ? a : b;
    for (unsigned i = 1; i <= count; ++i) {
        c -= safe_div(Tf(fa + (B + A * (c - b)) * (c - a)),
                      Tf(B + A * (2 * c - a - b)), Tf(1 + c - a));
    }

    return outside(c, a, b) ? secant_interpolate(a, b, fa, fb) : c;
}

template <typename Tx, typename Tf>
Tx cubic_interpolation(const Tx& a,
                       const Tx& b,
                       const Tx& d,
                       const Tx& e,
                       const Tf& fa,
                       const Tf& fb,
                       const Tf& fd,
                       const Tf& fe) {
    Tx q11 = (d - e) * fd / (fe - fd);
    Tx q21 = (b - d) * fb / (fd - fb);
    Tx q31 = (a - b) * fa / (fb - fa);
    Tx d21 = (b - d) * fd / (fd - fb);
    Tx d31 = (a - b) * fb / (fb - fa);

    Tx q22 = (d21 - q11) * fb / (fe - fb);
    Tx q32 = (d31 - q21) * fa / (fd - fa);
    Tx d32 = (d31 - q21) * fd / (fd - fa);
    Tx q33 = (d32 - q22) * fa / (fe - fa);

    Tx c = q31 + q32 + q33 + a;

    return outside(c, a, b) ? quadratic_interpolation(a, b, d, fa, fb, fd, 3)
                            : c;
}
}  // namespace detail

/**
 * @brief Find root of a function using TOMS748. A simplified version adapted
 * from Boost library, ignoring most sanity checks.
 *
 * @tparam T input value type
 * @param func
 * @param a left staring point
 * @param b right starting point
 */
template <typename Func, typename Tx>
Tx find_root(const Func& func, const Tx& ax, const Tx& bx) {
    using Tf = typename std::result_of<Func(Tx)>::type;
    // control parameters

    const unsigned max_iter = 50;
    unsigned count = max_iter;
    static const Tf mu = 0.5f;
    auto tol = [](const Tx& x, const Tx& y) {
        return util::abs(x - y) <=
               std::min(util::abs(x), util::abs(y)) *
                   std::sqrt(std::numeric_limits<Tf>::epsilon());
    };

    if (ax >= bx) { throw std::domain_error("Given interval do not exist."); }

    Tx a = ax, b = bx;
    Tf fa = func(a), fb = func(b);

    if (sgn(fa) * sgn(fb) > 0) {
        throw std::domain_error("Given interval may not bracket a root.");
    }

    if (tol(a, b) || (fa == 0) || (fb == 0)) { return fb == 0 ? b : a; }

    Tx d = 0;
    Tx e = 1e5F;
    Tf fe = 1e5F;
    Tf fd = 1e5F;

    Tx c = detail::secant_interpolate(a, b, fa, fb);
    detail::bracket(func, a, b, c, fa, fb, d, fd);
    --count;
    if (count != 0 && (fa != 0) && !(tol(a, b))) {
        c = detail::quadratic_interpolation(a, b, d, fa, fb, fd, 2);
        e = d;
        fe = fd;
        detail::bracket(func, a, b, c, fa, fb, d, fd);
        --count;
    }

    Tf min_diff = std::numeric_limits<Tf>::min() * 32;
    Tx u;
    Tf fu;
    while (count != 0 && (fa != 0) && !tol(a, b)) {
        Tx a0 = a;
        Tx b0 = b;

        bool prof = (util::abs(fa - fb) < min_diff) ||
                    (util::abs(fa - fd) < min_diff) ||
                    (util::abs(fa - fe) < min_diff) ||
                    (util::abs(fb - fd) < min_diff) ||
                    (util::abs(fb - fe) < min_diff) ||
                    (util::abs(fd - fe) < min_diff);
        c = prof ? detail::quadratic_interpolation(a, b, d, fa, fb, fd, 2)
                 : detail::cubic_interpolation(a, b, d, e, fa, fb, fd, fe);

        e = d;
        fe = fd;
        detail::bracket(func, a, b, c, fa, fb, d, fd);
        if ((0 == --count) || (fa == 0) || tol(a, b)) { break; }

        prof = (util::abs(fa - fb) < min_diff) ||
               (util::abs(fa - fd) < min_diff) ||
               (util::abs(fa - fe) < min_diff) ||
               (util::abs(fb - fd) < min_diff) ||
               (util::abs(fb - fe) < min_diff) ||
               (util::abs(fd - fe) < min_diff);
        c = prof ? detail::quadratic_interpolation(a, b, d, fa, fb, fd, 3)
                 : detail::cubic_interpolation(a, b, d, e, fa, fb, fd, fe);

        detail::bracket(func, a, b, c, fa, fb, d, fd);
        if ((0 == --count) || (fa == 0) || tol(a, b)) { break; }

        if (util::abs(fa) < util::abs(fb)) {
            u = a;
            fu = fa;
        } else {
            u = b;
            fu = fb;
        }
        c = u - 2 * (fu / (fb - fa)) * (b - a);
        if (util::abs(c - u) > (b - a) / 2) { c = a + (b - a) / 2; }

        e = d;
        fe = fd;
        detail::bracket(func, a, b, c, fa, fb, d, fd);
        if ((0 == --count) || (fa == 0) || tol(a, b)) { break; }

        if ((b - a) < mu * (b0 - a0)) { continue; }

        e = d;
        fe = fd;
        detail::bracket(func, a, b, Tx(a + (b - a) / 2), fa, fb, d, fd);
        --count;
    }

#ifdef _DEBUG
    std::cout << "Iterate " << max_iter - count << " times.\n";
#endif
    return fb == 0 ? b : a;
}
}  // namespace util
