#include <array>
#include <type_traits>
#include <utility>

/**
 * @brief Common functionality across all Vec type regardless of internal value
 * type and dimensionality.
 *
 * @tparam D Dimension of the underlying space.
 * @tparam T Type of coordinate.
 */
template <unsigned D, typename T>
class Vec_base;

/**
 * @brief A vector in Eucleadian space.
 *
 * @tparam D Dimension of the underlying space.
 * @tparam T Type of coordinate.
 */
template <unsigned D, typename T>
class Vec;

template <unsigned D, typename T>
class Vec_base {
    static_assert(D > 0, "D cannot be 0.");

   public:
    static constexpr unsigned dim = D;
    using value_type = T;
    using vec_type = Vec<dim, value_type>;

   protected:
    std::array<value_type, dim> coord;

   public:
    Vec_base() = default;

    template <typename... Ts,
              typename = typename std::enable_if<sizeof...(Ts) == dim>::type>
    Vec_base(Ts... vs) noexcept : coord{(value_type)vs...} {}

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
        for (unsigned i = 0; i < dim; ++i) {
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
        for (unsigned i = 0; i < dim; ++i) {
            if (coord[i] != rhs[i]) { return false; }
        }
        return true;
    }
    bool operator!=(const vec_type& rhs) const noexcept {
        return !operator==(rhs);
    }

    // arithmetic operations

    vec_type& operator+=(const vec_type& rhs) noexcept {
        for (unsigned i = 0; i < dim; ++i) { coord[i] += rhs[i]; }
        return *reinterpret_cast<vec_type*>(this);
    }
    vec_type& operator-=(const vec_type& rhs) noexcept {
        for (unsigned i = 0; i < dim; ++i) { coord[i] -= rhs[i]; }
        return *reinterpret_cast<vec_type*>(this);
    }

    vec_type& operator*=(const T& scalar) noexcept {
        for (unsigned i = 0; i < dim; ++i) { coord[i] *= scalar; }
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
};

template <unsigned D, typename T>
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
