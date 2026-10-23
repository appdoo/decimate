// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <cmath>
#include <cstddef>
#include <type_traits>

namespace decimate::math {

/// 4D vector (homogeneous coordinates, RGBA). Four contiguous T, standard layout, GPU-uploadable as-is.
///
/// Integer instantiations (Vec4i/Vec4u) support only the arithmetic and comparison members;
/// length(), normalize() and operator/ truncate and are not meaningful for them.
template <typename T>
struct Vec4 {
    static_assert(std::is_arithmetic<T>::value, "Vec4 requires arithmetic type");

    // Dual access: name the components (.x/.y/.z/.w) or index them (data[i]) —
    // both views alias the same four contiguous T (standard-layout; see asserts
    // below). Matches the Vec2/Vec3 union layout for family consistency.
    union {
        T data[4];
        struct {
            T x, y, z, w;
        };
    };

    /// Zero vector.
    constexpr Vec4() noexcept
        : x(0), y(0), z(0), w(0) {
    }
    /// Vector from components.
    constexpr Vec4(T x, T y, T z, T w) noexcept
        : x(x), y(y), z(z), w(w) {
    }
    /// All four components set to @p scalar.
    constexpr explicit Vec4(T scalar) noexcept
        : x(scalar), y(scalar), z(scalar), w(scalar) {
    }

    /// Cross-precision conversion (e.g. Vec4<double> → Vec4f); explicit
    /// so a change of element type is never silent.
    template <typename U>
    constexpr explicit Vec4(const Vec4<U> &o) noexcept
        : x(static_cast<T>(o.x)), y(static_cast<T>(o.y)),
          z(static_cast<T>(o.z)), w(static_cast<T>(o.w)) {
    }

    constexpr Vec4(const Vec4 &) noexcept = default;
    constexpr Vec4(Vec4 &&) noexcept = default;
    constexpr Vec4 &operator=(const Vec4 &) noexcept = default;
    constexpr Vec4 &operator=(Vec4 &&) noexcept = default;

    /// Component by index 0..3 via the union array; no bounds check. NOT constexpr:
    /// reading the inactive union member is ill-formed in a constant expression.
    T &operator[](size_t i) noexcept {
        return data[i];
    }
    const T &operator[](size_t i) const noexcept {
        return data[i];
    }

    /// Implicit pointer to the 4 contiguous components, for GPU upload / raw interop
    /// (the `data` member is equally accessible directly).
    operator T *() noexcept {
        return data;
    }
    operator const T *() const noexcept {
        return data;
    }

    /// Componentwise sum.
    constexpr Vec4 operator+(const Vec4 &v) const noexcept {
        return { x + v.x, y + v.y, z + v.z, w + v.w };
    }

    /// Componentwise difference.
    constexpr Vec4 operator-(const Vec4 &v) const noexcept {
        return { x - v.x, y - v.y, z - v.z, w - v.w };
    }

    /// Scale by @p scalar.
    constexpr Vec4 operator*(T scalar) const noexcept {
        return { x * scalar, y * scalar, z * scalar, w * scalar };
    }

    /// Divide by @p scalar; no zero check (float gives inf/NaN, integers are UB).
    constexpr Vec4 operator/(T scalar) const noexcept {
        return { x / scalar, y / scalar, z / scalar, w / scalar };
    }

    /// Negation.
    constexpr Vec4 operator-() const noexcept {
        return { -x, -y, -z, -w };
    }

    /// Componentwise (Hadamard) product, NOT a dot product; use dot() for that.
    constexpr Vec4 operator*(const Vec4 &v) const noexcept {
        return { x * v.x, y * v.y, z * v.z, w * v.w };
    }

    /// In-place componentwise sum.
    constexpr Vec4 &operator+=(const Vec4 &v) noexcept {
        x += v.x;
        y += v.y;
        z += v.z;
        w += v.w;
        return *this;
    }

    /// In-place componentwise difference.
    constexpr Vec4 &operator-=(const Vec4 &v) noexcept {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        w -= v.w;
        return *this;
    }

    /// In-place scale.
    constexpr Vec4 &operator*=(T scalar) noexcept {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        w *= scalar;
        return *this;
    }

    /// In-place divide; no zero check.
    constexpr Vec4 &operator/=(T scalar) noexcept {
        x /= scalar;
        y /= scalar;
        z /= scalar;
        w /= scalar;
        return *this;
    }

    /// In-place componentwise product.
    constexpr Vec4 &operator*=(const Vec4 &v) noexcept {
        x *= v.x;
        y *= v.y;
        z *= v.z;
        w *= v.w;
        return *this;
    }

    /// Exact componentwise comparison (no tolerance).
    constexpr bool operator==(const Vec4 &v) const noexcept {
        return x == v.x && y == v.y && z == v.z && w == v.w;
    }

    constexpr bool operator!=(const Vec4 &v) const noexcept {
        return !(*this == v);
    }

    /// Dot product.
    constexpr T dot(const Vec4 &v) const noexcept {
        return x * v.x + y * v.y + z * v.z + w * v.w;
    }

    /// Euclidean length.
    T length() const noexcept {
        return std::sqrt(x * x + y * y + z * z + w * w);
    }

    /// Squared length; avoids the sqrt.
    constexpr T lengthSquared() const noexcept {
        return x * x + y * y + z * z + w * w;
    }

    /// Euclidean distance to @p v.
    T distance(const Vec4 &v) const noexcept {
        return (*this - v).length();
    }

    /// Squared distance to @p v.
    constexpr T distanceSquared(const Vec4 &v) const noexcept {
        return (*this - v).lengthSquared();
    }

    /// Unit-length copy; the zero vector stays zero.
    Vec4 normalized() const noexcept {
        T len = length();
        return len > T(0) ? (*this / len) : Vec4 {};
    }

    /// Normalizes in place; the zero vector is left unchanged.
    void normalize() noexcept {
        T len = length();
        if (len > T(0)) {
            *this /= len;
        }
    }

    /// Linear interpolation toward @p v; t is not clamped, so t outside [0, 1] extrapolates.
    constexpr Vec4 lerp(const Vec4 &v, T t) const noexcept {
        return *this + (v - *this) * t;
    }

    /// Each component limited to [minVal, maxVal]; requires minVal <= maxVal.
    constexpr Vec4 clamp(T minVal, T maxVal) const noexcept {
        return {
            x < minVal ? minVal : (x > maxVal ? maxVal : x),
            y < minVal ? minVal : (y > maxVal ? maxVal : y),
            z < minVal ? minVal : (z > maxVal ? maxVal : z),
            w < minVal ? minVal : (w > maxVal ? maxVal : w)
        };
    }

    /// (0, 0, 0, 0).
    static constexpr Vec4 zero() noexcept {
        return { 0, 0, 0, 0 };
    }
    /// (1, 1, 1, 1).
    static constexpr Vec4 one() noexcept {
        return { 1, 1, 1, 1 };
    }
    /// (1, 0, 0, 0).
    static constexpr Vec4 unitX() noexcept {
        return { 1, 0, 0, 0 };
    }
    /// (0, 1, 0, 0).
    static constexpr Vec4 unitY() noexcept {
        return { 0, 1, 0, 0 };
    }
    /// (0, 0, 1, 0).
    static constexpr Vec4 unitZ() noexcept {
        return { 0, 0, 1, 0 };
    }
    /// (0, 0, 0, 1).
    static constexpr Vec4 unitW() noexcept {
        return { 0, 0, 0, 1 };
    }

    /// Componentwise minimum of @p a and @p b.
    static constexpr Vec4 min(const Vec4 &a, const Vec4 &b) noexcept {
        return {
            a.x < b.x ? a.x : b.x,
            a.y < b.y ? a.y : b.y,
            a.z < b.z ? a.z : b.z,
            a.w < b.w ? a.w : b.w
        };
    }

    /// Componentwise maximum of @p a and @p b.
    static constexpr Vec4 max(const Vec4 &a, const Vec4 &b) noexcept {
        return {
            a.x > b.x ? a.x : b.x,
            a.y > b.y ? a.y : b.y,
            a.z > b.z ? a.z : b.z,
            a.w > b.w ? a.w : b.w
        };
    }
};


/// Scalar times vector; commutes with Vec4::operator*(T).
template <typename T>
constexpr Vec4<T> operator*(T scalar, const Vec4<T> &v) noexcept {
    return v * scalar;
}


using Vec4f = Vec4<float>;
using Vec4d = Vec4<double>;
using Vec4i = Vec4<int>;
using Vec4u = Vec4<unsigned int>;


// Uploaded to the GPU as-is, so layout is part of the contract.
static_assert(std::is_standard_layout<Vec4f>::value, "Vec4f must be standard layout");
static_assert(sizeof(Vec4f) == sizeof(float) * 4, "Vec4f must be tightly packed");

} // namespace decimate::math
