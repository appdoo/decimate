// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <cmath>
#include <type_traits>

namespace decimate::math {

/// 2D vector, mainly for texture coordinates (UV). Two contiguous T, standard layout.
///
/// Integer instantiations (Vec2i/Vec2u) support only the arithmetic and comparison members;
/// length(), normalize() and operator/ truncate and are not meaningful for them.
template <typename T>
struct Vec2 {
    static_assert(std::is_arithmetic<T>::value, "Vec2 requires arithmetic type");

    // Dual access: name the components (.x/.y) or index them (data[i]) — both
    // views alias the same two contiguous T (standard-layout; see asserts
    // below). Matches the Vec3/Vec4 union layout for family consistency.
    // The anonymous struct inside the union is a compiler extension (clang/gcc/MSVC all accept it).
    union {
        T data[2];
        struct {
            T x, y;
        };
    };

    /// Zero vector.
    constexpr Vec2() noexcept
        : x(0), y(0) {
    }
    /// Vector from components.
    constexpr Vec2(T x, T y) noexcept
        : x(x), y(y) {
    }
    /// Both components set to @p scalar.
    constexpr explicit Vec2(T scalar) noexcept
        : x(scalar), y(scalar) {
    }

    /// Cross-precision conversion (e.g. Vec2<double> → Vec2f); explicit
    /// so a change of element type is never silent.
    template <typename U>
    constexpr explicit Vec2(const Vec2<U> &o) noexcept
        : x(static_cast<T>(o.x)), y(static_cast<T>(o.y)) {
    }

    constexpr Vec2(const Vec2 &) noexcept = default;
    constexpr Vec2(Vec2 &&) noexcept = default;
    constexpr Vec2 &operator=(const Vec2 &) noexcept = default;
    constexpr Vec2 &operator=(Vec2 &&) noexcept = default;

    /// Component by index 0..1 via the union array; no bounds check. NOT constexpr:
    /// reading the inactive union member is ill-formed in a constant expression.
    T &operator[](size_t i) noexcept {
        return data[i];
    }
    const T &operator[](size_t i) const noexcept {
        return data[i];
    }

    /// Implicit pointer to the 2 contiguous components, for GPU upload / raw interop
    /// (the `data` member is equally accessible directly).
    operator T *() noexcept {
        return data;
    }
    operator const T *() const noexcept {
        return data;
    }

    /// Componentwise sum.
    constexpr Vec2 operator+(const Vec2 &v) const noexcept {
        return { x + v.x, y + v.y };
    }

    /// Componentwise difference.
    constexpr Vec2 operator-(const Vec2 &v) const noexcept {
        return { x - v.x, y - v.y };
    }

    /// Scale by @p scalar.
    constexpr Vec2 operator*(T scalar) const noexcept {
        return { x * scalar, y * scalar };
    }

    /// Divide by @p scalar; no zero check (float gives inf/NaN, integers are UB).
    constexpr Vec2 operator/(T scalar) const noexcept {
        return { x / scalar, y / scalar };
    }

    /// Negation.
    constexpr Vec2 operator-() const noexcept {
        return { -x, -y };
    }

    /// In-place componentwise sum.
    constexpr Vec2 &operator+=(const Vec2 &v) noexcept {
        x += v.x;
        y += v.y;
        return *this;
    }

    /// In-place componentwise difference.
    constexpr Vec2 &operator-=(const Vec2 &v) noexcept {
        x -= v.x;
        y -= v.y;
        return *this;
    }

    /// In-place scale.
    constexpr Vec2 &operator*=(T scalar) noexcept {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    /// In-place divide; no zero check.
    constexpr Vec2 &operator/=(T scalar) noexcept {
        x /= scalar;
        y /= scalar;
        return *this;
    }

    /// Exact componentwise comparison (no tolerance).
    constexpr bool operator==(const Vec2 &v) const noexcept {
        return x == v.x && y == v.y;
    }

    constexpr bool operator!=(const Vec2 &v) const noexcept {
        return !(*this == v);
    }

    /// Dot product.
    constexpr T dot(const Vec2 &v) const noexcept {
        return x * v.x + y * v.y;
    }

    /// Euclidean length.
    T length() const noexcept {
        return std::sqrt(x * x + y * y);
    }

    /// Squared length; avoids the sqrt.
    constexpr T lengthSquared() const noexcept {
        return x * x + y * y;
    }

    /// Unit-length copy; the zero vector stays zero.
    Vec2 normalized() const noexcept {
        T len = length();
        return len > T(0) ? (*this / len) : Vec2 {};
    }

    /// Normalizes in place; the zero vector is left unchanged.
    void normalize() noexcept {
        T len = length();
        if (len > T(0)) {
            *this /= len;
        }
    }

    /// (0, 0).
    static constexpr Vec2 zero() noexcept {
        return { 0, 0 };
    }
    /// (1, 1).
    static constexpr Vec2 one() noexcept {
        return { 1, 1 };
    }
    /// (1, 0).
    static constexpr Vec2 unitX() noexcept {
        return { 1, 0 };
    }
    /// (0, 1).
    static constexpr Vec2 unitY() noexcept {
        return { 0, 1 };
    }
};


/// Scalar times vector; commutes with Vec2::operator*(T).
template <typename T>
constexpr Vec2<T> operator*(T scalar, const Vec2<T> &v) noexcept {
    return v * scalar;
}


using Vec2f = Vec2<float>;
using Vec2d = Vec2<double>;
using Vec2i = Vec2<int>;
using Vec2u = Vec2<unsigned int>;


// Uploaded to the GPU as-is, so layout is part of the contract.
static_assert(std::is_standard_layout<Vec2f>::value, "Vec2f must be standard layout");
static_assert(sizeof(Vec2f) == sizeof(float) * 2, "Vec2f must be tightly packed");

} // namespace decimate::math
