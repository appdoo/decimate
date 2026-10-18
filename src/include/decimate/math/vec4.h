// This file is part of the decimate library; see decimate.h for version/license details.
#pragma once

#include <cmath>
#include <type_traits>

namespace decimate::math {

/**
 * Vec4 - 4D vector template class
 *
 * Designed for high performance with:
 * - Contiguous memory layout for GPU upload (castable to T*)
 * - Move semantics and constexpr where possible
 * - Zero-overhead abstractions
 * - SIMD-friendly alignment (16 bytes for float)
 */
template<typename T>
struct Vec4 {
    static_assert(std::is_arithmetic<T>::value, "Vec4 requires arithmetic type");

    // Dual access: name the components (.x/.y/.z/.w) or index them (data[i]) —
    // both views alias the same four contiguous T (standard-layout; see asserts
    // below). Matches the Vec2/Vec3 union layout for family consistency.
    union {
        T data[4];
        struct { T x, y, z, w; };
    };

    // Constructors — the default zero-initializes.
    constexpr Vec4() noexcept : x(0), y(0), z(0), w(0) {}
    constexpr Vec4(T x, T y, T z, T w) noexcept : x(x), y(y), z(z), w(w) {}
    constexpr explicit Vec4(T scalar) noexcept : x(scalar), y(scalar), z(scalar), w(scalar) {}

    // Cross-precision conversion (e.g. Vec4<double> → Vec4<float>); explicit
    // so a change of element type is never silent.
    template<typename U>
    constexpr explicit Vec4(const Vec4<U>& o) noexcept
        : x(static_cast<T>(o.x)), y(static_cast<T>(o.y)),
          z(static_cast<T>(o.z)), w(static_cast<T>(o.w)) {}

    // Copy and move constructors (trivial for arithmetic T)
    constexpr Vec4(const Vec4&) noexcept = default;
    constexpr Vec4(Vec4&&) noexcept = default;
    constexpr Vec4& operator=(const Vec4&) noexcept = default;
    constexpr Vec4& operator=(Vec4&&) noexcept = default;

    // Array-like access via the union array. NOT constexpr: reading the
    // inactive union member is ill-formed in a constant expression.
    T& operator[](size_t i) noexcept { return data[i]; }
    const T& operator[](size_t i) const noexcept { return data[i]; }

    // Direct pointer access for GPU upload / raw interop (the `data` member is
    // equally accessible directly).
    operator T*() noexcept { return data; }
    operator const T*() const noexcept { return data; }

    // Arithmetic operators
    constexpr Vec4 operator+(const Vec4& v) const noexcept {
        return {x + v.x, y + v.y, z + v.z, w + v.w};
    }

    constexpr Vec4 operator-(const Vec4& v) const noexcept {
        return {x - v.x, y - v.y, z - v.z, w - v.w};
    }

    constexpr Vec4 operator*(T scalar) const noexcept {
        return {x * scalar, y * scalar, z * scalar, w * scalar};
    }

    constexpr Vec4 operator/(T scalar) const noexcept {
        return {x / scalar, y / scalar, z / scalar, w / scalar};
    }

    constexpr Vec4 operator-() const noexcept {
        return {-x, -y, -z, -w};
    }

    // Component-wise multiplication
    constexpr Vec4 operator*(const Vec4& v) const noexcept {
        return {x * v.x, y * v.y, z * v.z, w * v.w};
    }

    // Compound assignment operators
    constexpr Vec4& operator+=(const Vec4& v) noexcept {
        x += v.x; y += v.y; z += v.z; w += v.w;
        return *this;
    }

    constexpr Vec4& operator-=(const Vec4& v) noexcept {
        x -= v.x; y -= v.y; z -= v.z; w -= v.w;
        return *this;
    }

    constexpr Vec4& operator*=(T scalar) noexcept {
        x *= scalar; y *= scalar; z *= scalar; w *= scalar;
        return *this;
    }

    constexpr Vec4& operator/=(T scalar) noexcept {
        x /= scalar; y /= scalar; z /= scalar; w /= scalar;
        return *this;
    }

    constexpr Vec4& operator*=(const Vec4& v) noexcept {
        x *= v.x; y *= v.y; z *= v.z; w *= v.w;
        return *this;
    }

    // Comparison operators
    constexpr bool operator==(const Vec4& v) const noexcept {
        return x == v.x && y == v.y && z == v.z && w == v.w;
    }

    constexpr bool operator!=(const Vec4& v) const noexcept {
        return !(*this == v);
    }

    // Vector operations
    constexpr T dot(const Vec4& v) const noexcept {
        return x * v.x + y * v.y + z * v.z + w * v.w;
    }

    T length() const noexcept {
        return std::sqrt(x * x + y * y + z * z + w * w);
    }

    constexpr T lengthSquared() const noexcept {
        return x * x + y * y + z * z + w * w;
    }

    T distance(const Vec4& v) const noexcept {
        return (*this - v).length();
    }

    constexpr T distanceSquared(const Vec4& v) const noexcept {
        return (*this - v).lengthSquared();
    }

    Vec4 normalized() const noexcept {
        T len = length();
        return len > T(0) ? (*this / len) : Vec4{};
    }

    void normalize() noexcept {
        T len = length();
        if (len > T(0)) {
            *this /= len;
        }
    }

    // Linear interpolation
    constexpr Vec4 lerp(const Vec4& v, T t) const noexcept {
        return *this + (v - *this) * t;
    }

    // Clamp each component between min and max
    constexpr Vec4 clamp(T minVal, T maxVal) const noexcept {
        return {
            x < minVal ? minVal : (x > maxVal ? maxVal : x),
            y < minVal ? minVal : (y > maxVal ? maxVal : y),
            z < minVal ? minVal : (z > maxVal ? maxVal : z),
            w < minVal ? minVal : (w > maxVal ? maxVal : w)
        };
    }

    // Static factory methods
    static constexpr Vec4 zero() noexcept { return {0, 0, 0, 0}; }
    static constexpr Vec4 one() noexcept { return {1, 1, 1, 1}; }
    static constexpr Vec4 unitX() noexcept { return {1, 0, 0, 0}; }
    static constexpr Vec4 unitY() noexcept { return {0, 1, 0, 0}; }
    static constexpr Vec4 unitZ() noexcept { return {0, 0, 1, 0}; }
    static constexpr Vec4 unitW() noexcept { return {0, 0, 0, 1}; }

    // Min/max component-wise
    static constexpr Vec4 min(const Vec4& a, const Vec4& b) noexcept {
        return {
            a.x < b.x ? a.x : b.x,
            a.y < b.y ? a.y : b.y,
            a.z < b.z ? a.z : b.z,
            a.w < b.w ? a.w : b.w
        };
    }

    static constexpr Vec4 max(const Vec4& a, const Vec4& b) noexcept {
        return {
            a.x > b.x ? a.x : b.x,
            a.y > b.y ? a.y : b.y,
            a.z > b.z ? a.z : b.z,
            a.w > b.w ? a.w : b.w
        };
    }
};

// Non-member operators for scalar * vector
template<typename T>
constexpr Vec4<T> operator*(T scalar, const Vec4<T>& v) noexcept {
    return v * scalar;
}

// Common typedefs
using Vec4f = Vec4<float>;
using Vec4d = Vec4<double>;
using Vec4i = Vec4<int>;
using Vec4u = Vec4<unsigned int>;

// Ensure standard layout for GPU compatibility
static_assert(std::is_standard_layout<Vec4f>::value, "Vec4f must be standard layout");
static_assert(sizeof(Vec4f) == sizeof(float) * 4, "Vec4f must be tightly packed");

} // namespace decimate::math
