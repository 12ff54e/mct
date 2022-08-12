#pragma once

#include <array>
#include <cmath>
#include <type_traits>
#include <utility>

/**
 * @brief Common functionality across all Vec type regardless of internal value
 * type and dimensionality.
 *
 * @tparam D Dimension of the underlying space.
 * @tparam T Type of coordinate.
 */
template <std::size_t D, typename T>
class Vec_base;

/**
 * @brief A vector in Eucleadian space.
 *
 * @tparam D Dimension of the underlying space.
 * @tparam T Type of coordinate.
 */
template <std::size_t D, typename T>
class Vec;

template <std::size_t D, typename T>
class Vec_base {
    static_assert(D > 0, "D cannot be 0.");

   public:
    static constexpr std::size_t dim = D;
    using value_type = T;
    using vec_type = Vec<dim, value_type>;

   protected:
    std::array<value_type, dim> coord;

   public:
    constexpr Vec_base() = default;

    template <typename... Ts,
              typename = typename std::enable_if<sizeof...(Ts) == dim>::type>
    constexpr Vec_base(Ts... vs) noexcept : coord{(value_type)vs...} {}

    /**
     * @brief Conversion constructor from supported type
     *
     * @tparam U
     * @param other another vec
     */
    template <typename U,
              typename = typename std::enable_if<
                  std::is_convertible<U, T>::value>::type>
    Vec_base(Vec_base<dim, U> other) noexcept {
        for (std::size_t i = 0; i < dim; ++i) {
            coord[i] = static_cast<value_type>(other[i]);
        }
    }

    // zero element

    static constexpr vec_type zero() noexcept {
        return Vec_base<dim, value_type>{};
    }

    // element access

    T& operator[](int i) noexcept { return coord[i]; }

    const T& operator[](int i) const noexcept { return coord[i]; }

    // comparison operations

    bool operator==(const vec_type& rhs) const noexcept {
        for (std::size_t i = 0; i < dim; ++i) {
            if (coord[i] != rhs[i]) { return false; }
        }
        return true;
    }
    bool operator!=(const vec_type& rhs) const noexcept {
        return !operator==(rhs);
    }

    // arithmetic operations

    vec_type& operator+=(const vec_type& rhs) noexcept {
        for (std::size_t i = 0; i < dim; ++i) { coord[i] += rhs[i]; }
        return *reinterpret_cast<vec_type*>(this);
    }
    vec_type& operator-=(const vec_type& rhs) noexcept {
        for (std::size_t i = 0; i < dim; ++i) { coord[i] -= rhs[i]; }
        return *reinterpret_cast<vec_type*>(this);
    }

    vec_type& operator*=(const T& scalar) noexcept {
        for (std::size_t i = 0; i < dim; ++i) { coord[i] *= scalar; }
        return *reinterpret_cast<vec_type*>(this);
    }

    vec_type& operator/=(const T& scalar) noexcept {
        for (std::size_t i = 0; i < dim; ++i) { coord[i] /= scalar; }
        return *reinterpret_cast<vec_type*>(this);
    }

    friend vec_type operator+(vec_type lhs, const vec_type& rhs) noexcept {
        lhs += rhs;
        return lhs;
    }
    friend vec_type operator-(vec_type lhs, const vec_type& rhs) noexcept {
        lhs -= rhs;
        return lhs;
    }

    friend vec_type operator*(vec_type vec, const T& scalar) noexcept {
        return vec *= scalar;
    }
    friend vec_type operator*(const T& scalar, vec_type vec) noexcept {
        return vec *= scalar;
    }

    friend vec_type operator/(vec_type vec, const T& scalar) noexcept {
        return vec /= scalar;
    }

    // conversion operator

    // convert to the underlying coordinate array
    constexpr operator std::array<value_type, dim>() const { return coord; }

    // convert to derived class
    constexpr operator vec_type() { return static_cast<vec_type&&>(*this); }

    // properties

    T __L2_norm_square() const {
        T norm{};
        for (auto& c : coord) { norm += c * c; }
        return norm;
    }

    T mag() const { return std::sqrt(__L2_norm_square()); }
};

/**
 * @brief Generic Vector type of any dimension, default to consisting double as
 * coordinateds
 *
 * @tparam D Dimension
 * @tparam T Underlying types of each dimension
 */
template <std::size_t D, typename T = double>
class Vec : public Vec_base<D, T> {};

template <typename T>
class Vec<2, T> : public Vec_base<2, T> {
   public:
    using Vec_base<2, T>::Vec_base;

    T& x() noexcept { return this->coord[0]; }
    T& y() noexcept { return this->coord[1]; }

    const T& x() const noexcept { return this->coord[0]; }
    const T& y() const noexcept { return this->coord[1]; }
};

template <typename T>
class Vec<3, T> : public Vec_base<3, T> {
   public:
    using Vec_base<3, T>::Vec_base;

    T& x() noexcept { return this->coord[0]; }
    T& y() noexcept { return this->coord[1]; }
    T& z() noexcept { return this->coord[2]; }

    const T& x() const noexcept { return this->coord[0]; }
    const T& y() const noexcept { return this->coord[1]; }
    const T& z() const noexcept { return this->coord[2]; }
};

template <typename T = double>
Vec<2, T> cross(const Vec<2, T>& vec) {
    return Vec<2, T>{-vec.y(), vec.x()};
}

template <typename T = double>
Vec<3, T> cross(const Vec<3, T>& v1, const Vec<3, T>& v2) {
    return Vec<3, T>{v1.y() * v2.z() - v1.z() * v2.y(),
                     v1.z() * v2.x() - v1.x() * v2.z(),
                     v1.x() * v2.y() - v1.y() * v2.x()};
}

template <std::size_t D, typename T = double>
T dot(const Vec<D, T>& v1, const Vec<D, T>& v2) {
    T p = 0;
    for (std::size_t i = 0; i < D; ++i) { p += v1[i] * v2[i]; }
    return p;
}
