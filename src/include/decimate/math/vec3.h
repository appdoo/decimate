// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <cmath>
#include <type_traits>

namespace decimate::math {

/// 3D vector (positions, normals, directions). Three contiguous T, standard layout, so an array of
/// Vec3 can be uploaded to the GPU or reinterpreted as a flat T array.
///
/// Integer instantiations (Vec3i/Vec3u) support only the arithmetic and comparison members;
/// length(), normalize() and operator/ truncate and are not meaningful for them.
template <typename T>
struct Vec3 {
    static_assert(std::is_arithmetic<T>::value, "Vec3 requires arithmetic type");

    // Dual access: name the components (.x/.y/.z) or index them (data[i]) —
    // both views alias the same three contiguous T (standard-layout; see the
    // asserts at the bottom of this file). The anonymous struct inside the
    // union is a near-universal compiler extension (glm/Unreal rely on it) and
    // compiles clean on the project's clang/gcc/Emscripten toolchains.
    union {
        T data[3];
        struct {
            T x, y, z;
        };
    };

    /// Zero vector; callers depend on the zero-init (e.g. a std::vector<Vec3> normal
    /// accumulator that is then +='d into).
    constexpr Vec3() noexcept
        : x(0), y(0), z(0) {
    }
    /// Vector from components.
    constexpr Vec3(T x, T y, T z) noexcept
        : x(x), y(y), z(z) {
    }
    /// All three components set to @p scalar.
    constexpr explicit Vec3(T scalar) noexcept
        : x(scalar), y(scalar), z(scalar) {
    }

    /// Cross-precision conversion (e.g. Vec3<double> → Vec3f); explicit
    /// so a change of element type is never silent.
    template <typename U>
    constexpr explicit Vec3(const Vec3<U> &o) noexcept
        : x(static_cast<T>(o.x)), y(static_cast<T>(o.y)), z(static_cast<T>(o.z)) {
    }

    constexpr Vec3(const Vec3 &) noexcept = default;
    constexpr Vec3(Vec3 &&) noexcept = default;
    constexpr Vec3 &operator=(const Vec3 &) noexcept = default;
    constexpr Vec3 &operator=(Vec3 &&) noexcept = default;

    /// Component by index 0..2 via the union array; no bounds check. NOT constexpr:
    /// reading the inactive union member is ill-formed in a constant expression.
    T &operator[](size_t i) noexcept {
        return data[i];
    }
    const T &operator[](size_t i) const noexcept {
        return data[i];
    }

    /// Implicit pointer to the 3 contiguous components, for GPU upload / raw interop
    /// (the `data` member is equally accessible directly).
    operator T *() noexcept {
        return data;
    }
    operator const T *() const noexcept {
        return data;
    }

    /// Componentwise sum.
    constexpr Vec3 operator+(const Vec3 &v) const noexcept {
        return { x + v.x, y + v.y, z + v.z };
    }

    /// Componentwise difference.
    constexpr Vec3 operator-(const Vec3 &v) const noexcept {
        return { x - v.x, y - v.y, z - v.z };
    }

    /// Scale by @p scalar.
    constexpr Vec3 operator*(T scalar) const noexcept {
        return { x * scalar, y * scalar, z * scalar };
    }

    /// Divide by @p scalar; no zero check (float gives inf/NaN, integers are UB).
    constexpr Vec3 operator/(T scalar) const noexcept {
        return { x / scalar, y / scalar, z / scalar };
    }

    /// Negation.
    constexpr Vec3 operator-() const noexcept {
        return { -x, -y, -z };
    }

    /// In-place componentwise sum.
    constexpr Vec3 &operator+=(const Vec3 &v) noexcept {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    /// In-place componentwise difference.
    constexpr Vec3 &operator-=(const Vec3 &v) noexcept {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }

    /// In-place scale.
    constexpr Vec3 &operator*=(T scalar) noexcept {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    /// In-place divide; no zero check.
    constexpr Vec3 &operator/=(T scalar) noexcept {
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }

    /// Exact componentwise comparison (no tolerance).
    constexpr bool operator==(const Vec3 &v) const noexcept {
        return x == v.x && y == v.y && z == v.z;
    }

    constexpr bool operator!=(const Vec3 &v) const noexcept {
        return !(*this == v);
    }

    /// Dot product.
    constexpr T dot(const Vec3 &v) const noexcept {
        return x * v.x + y * v.y + z * v.z;
    }

    /// Cross product (right-handed): perpendicular to both, length |a||b|sin(angle).
    constexpr Vec3 cross(const Vec3 &v) const noexcept {
        return {
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        };
    }

    /// Componentwise (Hadamard) product, NOT a dot product; use dot() for that.
    constexpr Vec3 operator*(const Vec3 &v) const noexcept {
        return { x * v.x, y * v.y, z * v.z };
    }

    /// In-place componentwise product.
    constexpr Vec3 &operator*=(const Vec3 &v) noexcept {
        x *= v.x;
        y *= v.y;
        z *= v.z;
        return *this;
    }

    /// Euclidean length.
    T length() const noexcept {
        return std::sqrt(x * x + y * y + z * z);
    }

    /// Squared length; avoids the sqrt.
    constexpr T lengthSquared() const noexcept {
        return x * x + y * y + z * z;
    }

    /// Euclidean distance to @p v.
    T distance(const Vec3 &v) const noexcept {
        return (*this - v).length();
    }

    /// Squared distance to @p v.
    constexpr T distanceSquared(const Vec3 &v) const noexcept {
        return (*this - v).lengthSquared();
    }

    /// Unit-length copy; the zero vector stays zero.
    Vec3 normalized() const noexcept {
        T len = length();
        return len > T(0) ? (*this / len) : Vec3 {};
    }

    /// Normalizes in place; the zero vector is left unchanged.
    void normalize() noexcept {
        T len = length();
        if (len > T(0)) {
            *this /= len;
        }
    }

    /// Linear interpolation toward @p v; t is not clamped, so t outside [0, 1] extrapolates.
    constexpr Vec3 lerp(const Vec3 &v, T t) const noexcept {
        return *this + (v - *this) * t;
    }

    /// Mirror of this vector about the plane with unit @p normal (normal must be normalized).
    Vec3 reflect(const Vec3 &normal) const noexcept {
        return *this - normal * (T(2) * dot(normal));
    }

    /// Each component limited to [minVal, maxVal]; requires minVal <= maxVal.
    constexpr Vec3 clamp(T minVal, T maxVal) const noexcept {
        return {
            x < minVal ? minVal : (x > maxVal ? maxVal : x),
            y < minVal ? minVal : (y > maxVal ? maxVal : y),
            z < minVal ? minVal : (z > maxVal ? maxVal : z)
        };
    }

    /// (0, 0, 0).
    static constexpr Vec3 zero() noexcept {
        return { 0, 0, 0 };
    }
    /// (1, 1, 1).
    static constexpr Vec3 one() noexcept {
        return { 1, 1, 1 };
    }
    /// (1, 0, 0).
    static constexpr Vec3 unitX() noexcept {
        return { 1, 0, 0 };
    }
    /// (0, 1, 0).
    static constexpr Vec3 unitY() noexcept {
        return { 0, 1, 0 };
    }
    /// (0, 0, 1).
    static constexpr Vec3 unitZ() noexcept {
        return { 0, 0, 1 };
    }

    /// Componentwise minimum of @p a and @p b.
    static constexpr Vec3 min(const Vec3 &a, const Vec3 &b) noexcept {
        return {
            a.x < b.x ? a.x : b.x,
            a.y < b.y ? a.y : b.y,
            a.z < b.z ? a.z : b.z
        };
    }

    /// Componentwise maximum of @p a and @p b.
    static constexpr Vec3 max(const Vec3 &a, const Vec3 &b) noexcept {
        return {
            a.x > b.x ? a.x : b.x,
            a.y > b.y ? a.y : b.y,
            a.z > b.z ? a.z : b.z
        };
    }
};


/// Scalar times vector; commutes with Vec3::operator*(T).
template <typename T>
constexpr Vec3<T> operator*(T scalar, const Vec3<T> &v) noexcept {
    return v * scalar;
}


using Vec3f = Vec3<float>;
using Vec3d = Vec3<double>;
using Vec3i = Vec3<int>;
using Vec3u = Vec3<unsigned int>;


// Uploaded to the GPU as-is, so layout is part of the contract.
static_assert(std::is_standard_layout<Vec3f>::value, "Vec3f must be standard layout");
static_assert(sizeof(Vec3f) == sizeof(float) * 3, "Vec3f must be tightly packed");

} // namespace decimate::math
